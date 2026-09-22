from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

APP_DIR = Path(__file__).resolve().parent
DEFAULT_THEORY = APP_DIR / "mouse_20s_theoretical_masses.csv"
W_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def read_theory(uploaded_file=None) -> pd.DataFrame:
    if uploaded_file is None:
        theory = pd.read_csv(DEFAULT_THEORY)
    else:
        name = uploaded_file.name.lower()
        if name.endswith(".csv"):
            theory = pd.read_csv(uploaded_file)
        else:
            theory = pd.read_excel(uploaded_file)
    theory["Theoretical_Mass_Da"] = pd.to_numeric(
        theory["Theoretical_Mass_Da"], errors="coerce"
    )
    return theory.dropna(subset=["Theoretical_Mass_Da"]).reset_index(drop=True)


def docx_text_tokens(raw: bytes) -> list[str]:
    """Read visible Word text in document.xml, including ProSight text boxes/shapes."""
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    out = []
    for node in root.iter(f"{W_NS}t"):
        if node.text is not None:
            value = node.text.strip()
            if value:
                out.append(value)
    return out


def value_after(tokens: list[str], label: str, default: str = "") -> str:
    try:
        i = tokens.index(label)
        return tokens[i + 1] if i + 1 < len(tokens) else default
    except ValueError:
        return default


def parse_name_metadata(report_filename: str, data_filename: str = "") -> dict[str, str]:
    joined = f"{Path(report_filename).stem} {Path(data_filename).stem}".lower()

    age = ""
    if re.search(r"(?:^|[^a-z])(old|adult|aged)(?:[^a-z]|$)", joined):
        age = "Old/Adult"
    elif re.search(r"(?:^|[^a-z])(infant|young|juvenile|pup)(?:[^a-z]|$)", joined):
        age = "Infant/Young"

    organ = ""
    for key, label in {
        "liver": "Liver",
        "lung": "Lung",
        "brain": "Brain",
        "kidney": "Kidney",
        "kideny": "Kidney",
        "heart": "Heart",
        "spleen": "Spleen",
    }.items():
        if key in joined:
            organ = label
            break

    peak = ""
    # Peak number is intentionally taken from the DOCX report filename first.
    m = re.search(r"peak[\s_\-]*(\d+)", Path(report_filename).stem, flags=re.I)
    if m:
        peak = f"Peak {int(m.group(1))}"

    return {
        "Report_File": report_filename,
        "Peak": peak,
        "Age_Group": age,
        "Organ": organ,
    }


def parse_prosight_docx(raw: bytes, report_filename: str) -> tuple[pd.DataFrame, dict[str, str]]:
    tokens = docx_text_tokens(raw)

    data_file = value_after(tokens, "Data File Name")
    scan_number = value_after(tokens, "Scan Number")
    retention_time = value_after(tokens, "Retention Time")
    method_name = value_after(tokens, "Method Name")
    report_date = value_after(tokens, "Report Creation Date")

    meta = {
        **parse_name_metadata(report_filename, data_file),
        "Data_File_Name": data_file,
        "Retention_Time": retention_time,
        "Scan_Number": scan_number,
        "Method_Name": method_name,
        "Report_Creation_Date": report_date,
    }

    records = []
    # In ProSight Native reports, every deconvoluted result block begins with "Mass".
    i = 0
    while i < len(tokens):
        if tokens[i] != "Mass":
            i += 1
            continue
        if i + 1 >= len(tokens):
            break
        try:
            mass = float(tokens[i + 1].replace(",", ""))
        except ValueError:
            i += 1
            continue

        # Search only within this result block, stopping before the next "Mass" label.
        j = i + 2
        block = {}
        while j < len(tokens) and tokens[j] != "Mass":
            if tokens[j] in {
                "Charge State Range",
                "Sum Intensity",
                "Relative Intensity",
                "Top ID (Sequence)",
                "Top ID (PTMs)",
                "Mass Error (Da)",
                "Total # of IDs",
            }:
                key = tokens[j]
                # Empty ProSight fields may be followed immediately by the next label.
                if j + 1 < len(tokens) and tokens[j + 1] not in {
                    "Charge State Range",
                    "Sum Intensity",
                    "Relative Intensity",
                    "Top ID (Sequence)",
                    "Top ID (PTMs)",
                    "Mass Error (Da)",
                    "Total # of IDs",
                    "Mass",
                }:
                    block[key] = tokens[j + 1]
                else:
                    block[key] = ""
            j += 1

        sum_intensity = pd.to_numeric(block.get("Sum Intensity", ""), errors="coerce")
        rel_text = block.get("Relative Intensity", "")
        rel_num = pd.to_numeric(str(rel_text).replace("%", ""), errors="coerce")
        records.append(
            {
                **meta,
                "Detected_Mass_Da": mass,
                "Charge_State_Range": block.get("Charge State Range", ""),
                "Sum_Intensity": sum_intensity,
                "Relative_Intensity_pct": rel_num,
                "ProSight_Top_ID_Sequence": block.get("Top ID (Sequence)", ""),
                "ProSight_Top_ID_PTMs": block.get("Top ID (PTMs)", ""),
                "ProSight_Total_IDs": pd.to_numeric(block.get("Total # of IDs", ""), errors="coerce"),
                "Intensity_QC": "Kept (no intensity filtering)",
            }
        )
        i = max(j, i + 1)

    return pd.DataFrame(records), meta


def match_masses(df: pd.DataFrame, theory: pd.DataFrame, tolerance: float, unit: str) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    t = theory.copy().reset_index(drop=True)
    tm = t["Theoretical_Mass_Da"].to_numpy(dtype=float)
    rows = []

    for _, row in df.iterrows():
        detected = float(row["Detected_Mass_Da"])
        diffs = detected - tm
        idx = int(np.argmin(np.abs(diffs)))
        theo = float(tm[idx])
        delta_da = float(diffs[idx])
        delta_ppm = (delta_da / theo) * 1_000_000.0
        is_match = abs(delta_da) <= tolerance if unit == "Da" else abs(delta_ppm) <= tolerance
        hit = t.iloc[idx]

        rec = row.to_dict()
        rec.update(
            {
                "Status": "Matched / identified" if is_match else "Unmatched / possible new mass",
                "Nearest_Subunit": hit.get("Subunit", ""),
                "Nearest_Accession": hit.get("Accession", ""),
                "Nearest_Modification": hit.get("Modification", ""),
                "Nearest_Theoretical_Mass_Da": theo,
                "Distance_to_Nearest_Da": delta_da,
                "Distance_to_Nearest_ppm": delta_ppm,
                "Subunit": hit.get("Subunit", "") if is_match else "",
                "Accession": hit.get("Accession", "") if is_match else "",
                "Modification": hit.get("Modification", "") if is_match else "",
                "Theoretical_Mass_Da": theo if is_match else np.nan,
                "Error_Da": delta_da if is_match else np.nan,
                "Error_ppm": delta_ppm if is_match else np.nan,
            }
        )
        rows.append(rec)
    return pd.DataFrame(rows)


def report_plot(result: pd.DataFrame):
    p = result.copy()
    p["Label"] = np.where(
        p["Status"].eq("Matched / identified"),
        p["Subunit"].fillna("") + " — " + p["Modification"].fillna(""),
        "Unmatched",
    )
    fig = px.scatter(
        p,
        x="Detected_Mass_Da",
        y="Sum_Intensity",
        symbol="Status",
        size="Relative_Intensity_pct",
        hover_data={
            "Detected_Mass_Da": ":.2f",
            "Sum_Intensity": ":.3g",
            "Relative_Intensity_pct": True,
            "Charge_State_Range": True,
            "Label": True,
            "Status": True,
        },
        title=f"{p['Peak'].iloc[0] or 'LC peak'} — deconvoluted masses",
        labels={"Detected_Mass_Da": "Detected mass (Da)", "Sum_Intensity": "Sum intensity"},
    )
    fig.update_traces(marker={"line": {"width": 1}})
    fig.update_layout(height=420, legend_title_text="Assignment")
    return fig


def to_excel_bytes(all_results: pd.DataFrame, theory: pd.DataFrame, reports: pd.DataFrame) -> bytes:
    output = io.BytesIO()
    matched = all_results[all_results["Status"] == "Matched / identified"]
    unmatched = all_results[all_results["Status"] != "Matched / identified"]
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        all_results.to_excel(writer, index=False, sheet_name="All masses")
        matched.to_excel(writer, index=False, sheet_name="Matched_identified")
        unmatched.to_excel(writer, index=False, sheet_name="Unmatched_new_masses")
        reports.to_excel(writer, index=False, sheet_name="Report summary")
        theory.to_excel(writer, index=False, sheet_name="Theory used")
    return output.getvalue()


def load_docx_sources(mode: str):
    items: list[tuple[str, bytes]] = []

    if mode == "Upload DOCX files":
        uploaded = st.file_uploader(
            "Select all ProSight Native DOCX reports",
            type=["docx"],
            accept_multiple_files=True,
            help="You can select many DOCX reports at once.",
        )
        for f in uploaded or []:
            items.append((f.name, f.getvalue()))

    elif mode == "Upload a ZIP folder":
        z = st.file_uploader(
            "Upload a ZIP containing your DOCX reports",
            type=["zip"],
            help="Useful when an experiment folder contains many Peak1, Peak2, Peak3... reports.",
        )
        if z is not None:
            try:
                with zipfile.ZipFile(io.BytesIO(z.getvalue())) as zf:
                    for name in sorted(zf.namelist()):
                        if name.lower().endswith(".docx") and not Path(name).name.startswith("~$"):
                            items.append((Path(name).name, zf.read(name)))
            except zipfile.BadZipFile:
                st.error("That file is not a valid ZIP archive.")

    else:
        folder = st.text_input(
            "Local folder path",
            placeholder=r"C:\Users\you\Documents\ProSight_reports",
            help="Works only when Streamlit is running on the same computer that contains this folder.",
        )
        recursive = st.checkbox("Include DOCX files in subfolders", value=True)
        if folder:
            p = Path(folder).expanduser()
            if not p.exists() or not p.is_dir():
                st.error("Folder not found on the computer running Streamlit.")
            else:
                paths = sorted(p.rglob("*.docx") if recursive else p.glob("*.docx"))
                for path in paths:
                    if not path.name.startswith("~$"):
                        try:
                            items.append((path.name, path.read_bytes()))
                        except OSError as exc:
                            st.warning(f"Could not read {path.name}: {exc}")
    return items


st.set_page_config(page_title="Proteasome ProSight Mass Matcher", page_icon="🧬", layout="wide")
st.title("🧬 Proteasome ProSight Mass Matcher")
st.caption(
    "Batch-read ProSight Native DOCX reports, track organ/age/LC peak/source RAW file, "
    "match deconvoluted masses to your theoretical 20S list, preserve unmatched masses, and export the results."
)

with st.sidebar:
    st.header("Matching settings")
    tolerance_unit = st.radio("Tolerance unit", ["Da", "ppm"], horizontal=True)
    default_tol = 5.0 if tolerance_unit == "Da" else 200.0
    tolerance = st.number_input(
        f"Maximum matching error ({tolerance_unit})",
        min_value=0.0,
        value=default_tol,
        step=0.5 if tolerance_unit == "Da" else 10.0,
    )
    st.caption("Unmatched masses are never discarded; the nearest theoretical candidate is retained only as context.")

st.subheader("1. Theoretical proteasome mass list")
use_default = st.checkbox("Use bundled Mouse 20S theoretical mass list", value=True)
theory_upload = None
if not use_default:
    theory_upload = st.file_uploader("Upload theoretical list", type=["csv", "xlsx", "xls"], key="theory")
    if theory_upload is None:
        st.stop()
theory = read_theory(theory_upload)
if "Theoretical_Mass_Da" not in theory.columns:
    st.error("The theoretical list must contain a 'Theoretical_Mass_Da' column.")
    st.stop()
with st.expander(f"View theoretical list ({len(theory)} masses)"):
    st.dataframe(theory, use_container_width=True, hide_index=True)

st.subheader("2. Add ProSight Native reports")
mode = st.radio(
    "How do you want to provide the reports?",
    ["Upload DOCX files", "Upload a ZIP folder", "Read a local folder"],
    horizontal=True,
)
if mode == "Read a local folder":
    st.info(
        "Local-folder mode is for running the app on your own computer. A Streamlit Community Cloud app cannot directly browse folders on your laptop."
    )
sources = load_docx_sources(mode)

if not sources:
    st.info("Add ProSight DOCX reports to start the batch analysis.")
    st.stop()

results = []
report_summaries = []
parse_errors = []
for filename, raw in sources:
    try:
        parsed, meta = parse_prosight_docx(raw, filename)
        if parsed.empty:
            parse_errors.append(f"{filename}: no deconvoluted mass blocks were found")
            continue
        matched = match_masses(parsed, theory, tolerance, tolerance_unit)
        results.append(matched)
        report_summaries.append(
            {
                **meta,
                "Detected_Masses": len(matched),
                "Matched_Identified": int((matched["Status"] == "Matched / identified").sum()),
                "Unmatched": int((matched["Status"] != "Matched / identified").sum()),
            }
        )
    except Exception as exc:
        parse_errors.append(f"{filename}: {exc}")

if parse_errors:
    with st.expander(f"Files needing attention ({len(parse_errors)})"):
        for msg in parse_errors:
            st.warning(msg)

if not results:
    st.error("No usable ProSight masses were extracted from the supplied DOCX reports.")
    st.stop()

combined = pd.concat(results, ignore_index=True)
reports_df = pd.DataFrame(report_summaries)

st.subheader("3. Batch overview")
a, b, c, d = st.columns(4)
a.metric("DOCX reports", len(reports_df))
b.metric("Detected masses", len(combined))
c.metric("Matched / identified", int((combined["Status"] == "Matched / identified").sum()))
d.metric("Unmatched", int((combined["Status"] != "Matched / identified").sum()))
st.dataframe(reports_df, use_container_width=True, hide_index=True)

st.subheader("4. Each LC peak report")
for n, result in enumerate(results):
    label_bits = [result["Age_Group"].iloc[0], result["Organ"].iloc[0], result["Peak"].iloc[0]]
    label = " · ".join(x for x in label_bits if x) or result["Report_File"].iloc[0]
    with st.expander(f"{label} — {result['Report_File'].iloc[0]}", expanded=(n == 0)):
        st.caption(
            f"RAW: {result['Data_File_Name'].iloc[0] or 'not found'}   |   "
            f"RT: {result['Retention_Time'].iloc[0] or 'not found'}"
        )
        left, right = st.columns([1.15, 1])
        with left:
            view_cols = [
                "Detected_Mass_Da", "Charge_State_Range", "Sum_Intensity", "Relative_Intensity_pct",
                "Status", "Subunit", "Modification", "Theoretical_Mass_Da", "Error_Da", "Error_ppm",
                "Nearest_Subunit", "Nearest_Modification", "Distance_to_Nearest_Da",
            ]
            st.dataframe(result[view_cols], use_container_width=True, hide_index=True)
        with right:
            st.plotly_chart(report_plot(result), use_container_width=True, key=f"plot_{n}")

st.subheader("5. Combined results")
tab_all, tab_matched, tab_unmatched = st.tabs(["All masses", "Matched / identified", "Unmatched / possible new masses"])
matched_only = combined[combined["Status"] == "Matched / identified"].copy()
unmatched_only = combined[combined["Status"] != "Matched / identified"].copy()
with tab_all:
    st.dataframe(combined, use_container_width=True, hide_index=True)
with tab_matched:
    st.dataframe(matched_only, use_container_width=True, hide_index=True)
with tab_unmatched:
    st.dataframe(unmatched_only, use_container_width=True, hide_index=True)
    st.caption("These masses remain in the dataset for later PTM/proteoform investigation; they are not force-assigned.")

st.subheader("6. Download")
excel = to_excel_bytes(combined, theory, reports_df)
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.download_button(
        "Download full Excel workbook",
        data=excel,
        file_name="proteasome_prosight_mass_analysis.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
with c2:
    st.download_button(
        "Download all masses CSV",
        data=combined.to_csv(index=False).encode("utf-8"),
        file_name="all_masses.csv",
        mime="text/csv",
    )
with c3:
    st.download_button(
        "Download matched CSV",
        data=matched_only.to_csv(index=False).encode("utf-8"),
        file_name="matched_identified_masses.csv",
        mime="text/csv",
    )
with c4:
    st.download_button(
        "Download unmatched CSV",
        data=unmatched_only.to_csv(index=False).encode("utf-8"),
        file_name="unmatched_possible_new_masses.csv",
        mime="text/csv",
    )
