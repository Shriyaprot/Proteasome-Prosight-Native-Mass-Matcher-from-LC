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
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

APP_DIR = Path(__file__).resolve().parent
ROOT_THEORY = APP_DIR / "mouse_20s_theoretical_masses.csv"
DATA_THEORY = APP_DIR / "data" / "mouse_20s_theoretical_masses.csv"
DEFAULT_THEORY = ROOT_THEORY if ROOT_THEORY.exists() else DATA_THEORY
W_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

# Putative neutral/adduct/PTM mass shifts used only for SECONDARY annotation.
# Direct theoretical matches are always attempted first.
COMMON_SHIFTS = pd.DataFrame(
    [
        ("Oxidation", 15.994915, "PTM"),
        ("Acetylation", 42.010565, "PTM"),
        ("Methylation", 14.015650, "PTM"),
        ("Dimethylation", 28.031300, "PTM"),
        ("Phosphorylation", 79.966331, "PTM"),
        ("Cysteinylation", 119.004099, "PTM"),
        ("Glutathionylation", 305.068156, "PTM"),
        ("Ammonia / ammonium-related", 17.026549, "Adduct"),
        ("Acetate", 59.013304, "Adduct"),
        ("Ammonium acetate", 77.039853, "Adduct"),
        ("Sodium replacing proton", 21.981943, "Adduct"),
        ("Potassium replacing proton", 37.955882, "Adduct"),
        ("Formate", 44.998201, "Adduct"),
        ("Acetonitrile", 41.026549, "Adduct"),
        ("DMSO", 78.013936, "Adduct"),
    ],
    columns=["Shift_Name", "Shift_Da", "Shift_Type"],
)

RESULT_COLUMNS = [
    "Detected_Mass_Da",
    "Charge_State_Range",
    "Sum_Intensity",
    "Relative_Intensity_pct",
    "Category",
    "Subunit",
    "Accession",
    "Modification",
    "Theoretical_Mass_Da",
    "Error_Da",
    "Error_ppm",
    "Possible_Shift_Type",
    "Possible_Shift",
    "Observed_Shift_Da",
    "Shift_Residual_Da",
    "Shift_Residual_ppm",
    "Reference_Subunit",
    "Reference_Modification",
    "Reference_Theoretical_Mass_Da",
]


def read_theory(uploaded_file=None) -> pd.DataFrame:
    if uploaded_file is None:
        theory = pd.read_csv(DEFAULT_THEORY)
    else:
        name = uploaded_file.name.lower()
        theory = pd.read_csv(uploaded_file) if name.endswith(".csv") else pd.read_excel(uploaded_file)

    # Keep only useful named fields; this also prevents empty Excel columns such as V1/W1/X1
    # or Unnamed:* columns from propagating into the output.
    aliases = {
        "Calculated Mass (Da)": "Theoretical_Mass_Da",
        "Calculated_Mass_Da": "Theoretical_Mass_Da",
        "Theoretical Mass (Da)": "Theoretical_Mass_Da",
    }
    theory = theory.rename(columns=aliases)
    wanted = [c for c in ["Subunit", "Accession", "Modification", "Theoretical_Mass_Da"] if c in theory.columns]
    if "Theoretical_Mass_Da" not in wanted:
        raise ValueError("The theoretical list needs a Theoretical_Mass_Da or Calculated Mass (Da) column.")
    theory = theory[wanted].copy()
    for c in ["Subunit", "Accession", "Modification"]:
        if c not in theory.columns:
            theory[c] = ""
        theory[c] = theory[c].ffill().fillna("")
    theory["Theoretical_Mass_Da"] = pd.to_numeric(theory["Theoretical_Mass_Da"], errors="coerce")
    return theory.dropna(subset=["Theoretical_Mass_Da"]).reset_index(drop=True)


def docx_text_tokens(raw: bytes) -> list[str]:
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
        "liver": "Liver", "lung": "Lung", "brain": "Brain", "kidney": "Kidney",
        "kideny": "Kidney", "heart": "Heart", "spleen": "Spleen",
    }.items():
        if key in joined:
            organ = label
            break

    peak = ""
    m = re.search(r"peak[\s_\-]*(\d+)", Path(report_filename).stem, flags=re.I)
    if m:
        peak = f"Peak {int(m.group(1))}"
    return {"Report_File": report_filename, "Peak": peak, "Age_Group": age, "Organ": organ}


def parse_charge_range(value: str) -> tuple[float, float, str]:
    text = str(value or "").strip().replace("–", "-").replace("—", "-")
    m = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", text)
    if not m:
        return np.nan, np.nan, text
    lo, hi = int(m.group(1)), int(m.group(2))
    return lo, hi, f"{lo} to {hi}"


def parse_prosight_docx(raw: bytes, report_filename: str) -> tuple[pd.DataFrame, dict[str, str]]:
    tokens = docx_text_tokens(raw)
    data_file = value_after(tokens, "Data File Name")
    meta = {
        **parse_name_metadata(report_filename, data_file),
        "Data_File_Name": data_file,
        "Retention_Time": value_after(tokens, "Retention Time"),
        "Scan_Number": value_after(tokens, "Scan Number"),
        "Method_Name": value_after(tokens, "Method Name"),
        "Report_Creation_Date": value_after(tokens, "Report Creation Date"),
    }

    labels = {
        "Charge State Range", "Sum Intensity", "Relative Intensity", "Top ID (Sequence)",
        "Top ID (PTMs)", "Mass Error (Da)", "Total # of IDs",
    }
    records = []
    i = 0
    while i < len(tokens):
        if tokens[i] != "Mass" or i + 1 >= len(tokens):
            i += 1
            continue
        try:
            mass = float(tokens[i + 1].replace(",", ""))
        except ValueError:
            i += 1
            continue

        j = i + 2
        block = {}
        while j < len(tokens) and tokens[j] != "Mass":
            if tokens[j] in labels:
                key = tokens[j]
                if j + 1 < len(tokens) and tokens[j + 1] not in labels | {"Mass"}:
                    block[key] = tokens[j + 1]
                else:
                    block[key] = ""
            j += 1

        lo, hi, safe_range = parse_charge_range(block.get("Charge State Range", ""))
        records.append({
            "Detected_Mass_Da": mass,
            "Charge_State_Min": lo,
            "Charge_State_Max": hi,
            "Charge_State_Range": safe_range,
            "Sum_Intensity": pd.to_numeric(block.get("Sum Intensity", ""), errors="coerce"),
            "Relative_Intensity_pct": pd.to_numeric(str(block.get("Relative Intensity", "")).replace("%", ""), errors="coerce"),
        })
        i = max(j, i + 1)

    return pd.DataFrame(records), meta


def ppm_error(observed: float, expected: float) -> float:
    return (observed - expected) / expected * 1_000_000.0


def match_masses(df: pd.DataFrame, theory: pd.DataFrame, direct_ppm: float, shift_ppm: float) -> pd.DataFrame:
    """Classify direct theoretical matches first, then putative PTM/adduct shifts.

    A mass outside the direct tolerance is NOT automatically called a PTM. It becomes
    Possible PTM/adduct only when its mass difference from a theoretical proteoform
    agrees with one of the configured shift masses. Otherwise it remains Unassigned.
    """
    if df.empty:
        return df.copy()

    t = theory.reset_index(drop=True)
    tm = t["Theoretical_Mass_Da"].to_numpy(dtype=float)
    shifts = COMMON_SHIFTS.reset_index(drop=True)
    rows = []

    for _, row in df.iterrows():
        detected = float(row["Detected_Mass_Da"])

        direct_ppms = (detected - tm) / tm * 1_000_000.0
        direct_idx = int(np.argmin(np.abs(direct_ppms)))
        direct_hit = t.iloc[direct_idx]
        direct_theo = float(tm[direct_idx])
        direct_da = detected - direct_theo
        direct_err_ppm = float(direct_ppms[direct_idx])

        rec = row.to_dict()
        rec.update({
            "Category": "Unassigned / investigate",
            "Subunit": "",
            "Accession": "",
            "Modification": "",
            "Theoretical_Mass_Da": np.nan,
            "Error_Da": np.nan,
            "Error_ppm": np.nan,
            "Possible_Shift_Type": "",
            "Possible_Shift": "",
            "Observed_Shift_Da": np.nan,
            "Shift_Residual_Da": np.nan,
            "Shift_Residual_ppm": np.nan,
            "Reference_Subunit": direct_hit.get("Subunit", ""),
            "Reference_Modification": direct_hit.get("Modification", ""),
            "Reference_Theoretical_Mass_Da": direct_theo,
        })

        if abs(direct_err_ppm) <= direct_ppm:
            rec.update({
                "Category": "Matched",
                "Subunit": direct_hit.get("Subunit", ""),
                "Accession": direct_hit.get("Accession", ""),
                "Modification": direct_hit.get("Modification", ""),
                "Theoretical_Mass_Da": direct_theo,
                "Error_Da": direct_da,
                "Error_ppm": direct_err_ppm,
                "Reference_Subunit": "",
                "Reference_Modification": "",
                "Reference_Theoretical_Mass_Da": np.nan,
            })
            rows.append(rec)
            continue

        # Search all theory x common shifts, allowing both + and - shift directions.
        best = None
        for ti, theo in enumerate(tm):
            observed_shift = detected - float(theo)
            for _, s in shifts.iterrows():
                shift_da = float(s["Shift_Da"])
                for direction in (1.0, -1.0):
                    expected = float(theo) + direction * shift_da
                    resid_da = detected - expected
                    resid_ppm = ppm_error(detected, expected)
                    score = abs(resid_ppm)
                    if best is None or score < best[0]:
                        best = (score, ti, s, direction, observed_shift, resid_da, resid_ppm, expected)

        if best is not None and best[0] <= shift_ppm:
            _, ti, s, direction, observed_shift, resid_da, resid_ppm, expected = best
            ref = t.iloc[ti]
            sign = "+" if direction > 0 else "−"
            rec.update({
                "Category": "Possible PTM" if s["Shift_Type"] == "PTM" else "Possible adduct",
                "Possible_Shift_Type": s["Shift_Type"],
                "Possible_Shift": f"{sign}{s['Shift_Name']} ({sign}{float(s['Shift_Da']):.4f} Da)",
                "Observed_Shift_Da": observed_shift,
                "Shift_Residual_Da": resid_da,
                "Shift_Residual_ppm": resid_ppm,
                "Reference_Subunit": ref.get("Subunit", ""),
                "Reference_Modification": ref.get("Modification", ""),
                "Reference_Theoretical_Mass_Da": float(tm[ti]),
            })

        rows.append(rec)
    return pd.DataFrame(rows)


def simple_result_view(result: pd.DataFrame) -> pd.DataFrame:
    return result[[c for c in RESULT_COLUMNS if c in result.columns]].copy()


def report_plot(result: pd.DataFrame):
    p = result.copy()
    p["Assignment"] = np.select(
        [
            p["Category"].eq("Matched"),
            p["Category"].eq("Possible PTM"),
            p["Category"].eq("Possible adduct"),
        ],
        [
            p["Subunit"].fillna("") + " — " + p["Modification"].fillna(""),
            "Possible PTM: " + p["Possible_Shift"].fillna(""),
            "Possible adduct: " + p["Possible_Shift"].fillna(""),
        ],
        default="Unassigned",
    )
    fig = px.scatter(
        p, x="Detected_Mass_Da", y="Sum_Intensity", symbol="Category",
        size="Relative_Intensity_pct",
        hover_data={
            "Detected_Mass_Da": ":.2f", "Sum_Intensity": ":.3g",
            "Relative_Intensity_pct": True, "Charge_State_Range": True,
            "Assignment": True, "Category": True,
        },
        title="Deconvoluted masses",
        labels={"Detected_Mass_Da": "Detected mass (Da)", "Sum_Intensity": "Sum intensity"},
    )
    fig.update_traces(marker={"line": {"width": 1}})
    fig.update_layout(height=420, legend_title_text="Classification")
    return fig


def style_sheet(ws, freeze: str | None = None):
    ws.sheet_view.showGridLines = False
    if freeze:
        ws.freeze_panes = freeze
    for cell in ws[1]:
        cell.font = Font(bold=True)
    widths = {}
    for row in ws.iter_rows():
        for cell in row:
            if cell.value is not None:
                widths[cell.column] = min(max(widths.get(cell.column, 0), len(str(cell.value)) + 2), 45)
    for idx, width in widths.items():
        ws.column_dimensions[get_column_letter(idx)].width = width


def append_result_table(ws, start_row: int, result: pd.DataFrame) -> int:
    view = simple_result_view(result)
    headers = list(view.columns)
    for col, h in enumerate(headers, 1):
        c = ws.cell(start_row, col, h)
        c.font = Font(bold=True)
        c.fill = PatternFill("solid", fgColor="D9EAF7")
        c.alignment = Alignment(wrap_text=True)
    for r_idx, values in enumerate(view.itertuples(index=False, name=None), start_row + 1):
        for c_idx, value in enumerate(values, 1):
            if pd.isna(value):
                value = None
            cell = ws.cell(r_idx, c_idx, value)
            if headers[c_idx - 1] == "Charge_State_Range":
                cell.number_format = "@"
    return start_row + len(view) + 1


def to_excel_bytes(results: list[pd.DataFrame], metas: list[dict[str, str]], theory: pd.DataFrame) -> bytes:
    """Create a simple grouped workbook: report metadata once, then mass rows."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Combined results"
    ws.sheet_view.showGridLines = False

    title_fill = PatternFill("solid", fgColor="1F4E78")
    title_font = Font(color="FFFFFF", bold=True, size=12)
    meta_fill = PatternFill("solid", fgColor="EAF2F8")

    row = 1
    for result, meta in zip(results, metas):
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=8)
        tc = ws.cell(row, 1, f"{meta.get('Report_File','')}  |  {meta.get('Peak','')}")
        tc.fill = title_fill
        tc.font = title_font
        row += 1

        meta_pairs = [
            ("Report file", meta.get("Report_File", "")),
            ("Data file", meta.get("Data_File_Name", "")),
            ("Age group", meta.get("Age_Group", "")),
            ("Organ", meta.get("Organ", "")),
            ("LC peak", meta.get("Peak", "")),
            ("Retention time", meta.get("Retention_Time", "")),
            ("Scan number", meta.get("Scan_Number", "")),
            ("Method", meta.get("Method_Name", "")),
            ("Report date", meta.get("Report_Creation_Date", "")),
        ]
        for label, value in meta_pairs:
            ws.cell(row, 1, label).font = Font(bold=True)
            ws.cell(row, 1).fill = meta_fill
            ws.cell(row, 2, value)
            row += 1
        row += 1
        row = append_result_table(ws, row, result)
        row += 2

    # Widths on combined sheet.
    for idx in range(1, ws.max_column + 1):
        max_len = 0
        for cell in ws[get_column_letter(idx)]:
            if cell.value is not None:
                max_len = max(max_len, min(len(str(cell.value)), 45))
        ws.column_dimensions[get_column_letter(idx)].width = max(12, min(max_len + 2, 45))

    # Category sheets are deliberately flat and minimal.
    combined = pd.concat(results, ignore_index=True)
    for sheet_name, category in [
        ("Matched", "Matched"),
        ("Possible PTM", "Possible PTM"),
        ("Possible adduct", "Possible adduct"),
        ("Unassigned", "Unassigned / investigate"),
    ]:
        sub = simple_result_view(combined[combined["Category"] == category])
        ws2 = wb.create_sheet(sheet_name)
        for c_idx, h in enumerate(sub.columns, 1):
            cell = ws2.cell(1, c_idx, h)
            cell.font = Font(bold=True)
            cell.fill = PatternFill("solid", fgColor="D9EAF7")
        for r_idx, vals in enumerate(sub.itertuples(index=False, name=None), 2):
            for c_idx, val in enumerate(vals, 1):
                if pd.isna(val):
                    val = None
                cell = ws2.cell(r_idx, c_idx, val)
                if sub.columns[c_idx - 1] == "Charge_State_Range":
                    cell.number_format = "@"
        style_sheet(ws2, freeze="A2")

    ws3 = wb.create_sheet("Theory used")
    clean_theory = theory[["Subunit", "Accession", "Modification", "Theoretical_Mass_Da"]]
    for c_idx, h in enumerate(clean_theory.columns, 1):
        cell = ws3.cell(1, c_idx, h)
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="D9EAD3")
    for r_idx, vals in enumerate(clean_theory.itertuples(index=False, name=None), 2):
        for c_idx, val in enumerate(vals, 1):
            ws3.cell(r_idx, c_idx, None if pd.isna(val) else val)
    style_sheet(ws3, freeze="A2")

    ws4 = wb.create_sheet("Shift reference")
    note = "Putative secondary annotations only. A mass is called Possible PTM/adduct only when the residual after applying the shift is within the configured shift tolerance."
    ws4.cell(1, 1, note)
    ws4.merge_cells(start_row=1, start_column=1, end_row=1, end_column=4)
    ws4.cell(1, 1).alignment = Alignment(wrap_text=True)
    for c_idx, h in enumerate(COMMON_SHIFTS.columns, 1):
        cell = ws4.cell(3, c_idx, h)
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="FFF2CC")
    for r_idx, vals in enumerate(COMMON_SHIFTS.itertuples(index=False, name=None), 4):
        for c_idx, val in enumerate(vals, 1):
            ws4.cell(r_idx, c_idx, val)
    style_sheet(ws4)

    output = io.BytesIO()
    wb.save(output)
    return output.getvalue()


def csv_bytes(df: pd.DataFrame) -> bytes:
    safe = simple_result_view(df)
    # UTF-8 with BOM opens cleanly in Excel; range text uses "to" so it cannot become a date.
    return safe.to_csv(index=False).encode("utf-8-sig")


def load_docx_sources(mode: str):
    items: list[tuple[str, bytes]] = []
    if mode == "Upload DOCX files":
        uploaded = st.file_uploader("Select all ProSight Native DOCX reports", type=["docx"], accept_multiple_files=True)
        for f in uploaded or []:
            items.append((f.name, f.getvalue()))
    elif mode == "Upload a ZIP folder":
        z = st.file_uploader("Upload a ZIP containing your DOCX reports", type=["zip"])
        if z is not None:
            try:
                with zipfile.ZipFile(io.BytesIO(z.getvalue())) as zf:
                    for name in sorted(zf.namelist()):
                        if name.lower().endswith(".docx") and not Path(name).name.startswith("~$"):
                            items.append((Path(name).name, zf.read(name)))
            except zipfile.BadZipFile:
                st.error("That file is not a valid ZIP archive.")
    else:
        folder = st.text_input("Local folder path", placeholder=r"C:\Users\you\Documents\ProSight_reports")
        recursive = st.checkbox("Include DOCX files in subfolders", value=True)
        if folder:
            p = Path(folder).expanduser()
            if not p.exists() or not p.is_dir():
                st.error("Folder not found on the computer running Streamlit.")
            else:
                paths = sorted(p.rglob("*.docx") if recursive else p.glob("*.docx"))
                for path in paths:
                    if not path.name.startswith("~$"):
                        items.append((path.name, path.read_bytes()))
    return items


st.set_page_config(page_title="Proteasome ProSight Mass Matcher", page_icon="🧬", layout="wide")
st.title("🧬 Proteasome ProSight Mass Matcher")
st.caption("Batch-read ProSight Native DOCX reports, match intact masses to your proteasome mass list, flag plausible PTM/adduct shifts, and retain unassigned masses.")

with st.sidebar:
    st.header("Matching settings")
    direct_ppm = st.number_input("Direct theoretical match tolerance (ppm)", min_value=0.1, value=5.0, step=0.5)
    shift_ppm = st.number_input("PTM/adduct residual tolerance (ppm)", min_value=0.1, value=5.0, step=0.5)
    st.caption("Direct matches are tested first. Outside that tolerance, PTM/adduct labels are only assigned when a configured mass shift also fits within the residual tolerance.")

st.subheader("1. Theoretical proteasome mass list")
use_default = st.checkbox("Use bundled Mouse 20S theoretical mass list", value=True)
theory_upload = None
if not use_default:
    theory_upload = st.file_uploader("Upload theoretical list", type=["csv", "xlsx", "xls"], key="theory")
    if theory_upload is None:
        st.stop()
try:
    theory = read_theory(theory_upload)
except Exception as exc:
    st.error(str(exc))
    st.stop()
with st.expander(f"View theoretical list ({len(theory)} masses)"):
    st.dataframe(theory, use_container_width=True, hide_index=True)

st.subheader("2. Add ProSight Native reports")
mode = st.radio("How do you want to provide the reports?", ["Upload DOCX files", "Upload a ZIP folder", "Read a local folder"], horizontal=True)
if mode == "Read a local folder":
    st.info("Local-folder mode works only when Streamlit is running on the same computer as the folder. Streamlit Community Cloud cannot browse a folder on your laptop.")
sources = load_docx_sources(mode)
if not sources:
    st.info("Add ProSight DOCX reports to start the analysis.")
    st.stop()

results: list[pd.DataFrame] = []
metas: list[dict[str, str]] = []
parse_errors = []
for filename, raw in sources:
    try:
        parsed, meta = parse_prosight_docx(raw, filename)
        if parsed.empty:
            parse_errors.append(f"{filename}: no deconvoluted mass blocks were found")
            continue
        result = match_masses(parsed, theory, direct_ppm, shift_ppm)
        results.append(result)
        metas.append(meta)
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

st.subheader("3. Batch overview")
counts = combined["Category"].value_counts()
a, b, c, d, e = st.columns(5)
a.metric("DOCX reports", len(results))
b.metric("Detected masses", len(combined))
c.metric("Matched", int(counts.get("Matched", 0)))
d.metric("Possible PTM/adduct", int(counts.get("Possible PTM", 0) + counts.get("Possible adduct", 0)))
e.metric("Unassigned", int(counts.get("Unassigned / investigate", 0)))

st.subheader("4. Each LC peak report")
for n, (result, meta) in enumerate(zip(results, metas)):
    label_bits = [meta.get("Age_Group", ""), meta.get("Organ", ""), meta.get("Peak", "")]
    label = " · ".join(x for x in label_bits if x) or meta.get("Report_File", "")
    with st.expander(f"{label} — {meta.get('Report_File','')}", expanded=(n == 0)):
        st.markdown(
            f"**Data file:** `{meta.get('Data_File_Name','') or 'not found'}`  \n"
            f"**Retention time:** {meta.get('Retention_Time','') or 'not found'}  |  "
            f"**Scans:** {meta.get('Scan_Number','') or 'not found'}"
        )
        left, right = st.columns([1.2, 1])
        with left:
            st.dataframe(simple_result_view(result), use_container_width=True, hide_index=True)
        with right:
            st.plotly_chart(report_plot(result), use_container_width=True, key=f"plot_{n}")

st.subheader("5. Combined mass results")
tab_all, tab_match, tab_shift, tab_unassigned = st.tabs(["All masses", "Matched", "Possible PTM/adduct", "Unassigned"])
matched = combined[combined["Category"] == "Matched"].copy()
possible = combined[combined["Category"].isin(["Possible PTM", "Possible adduct"])].copy()
unassigned = combined[combined["Category"] == "Unassigned / investigate"].copy()
with tab_all:
    st.dataframe(simple_result_view(combined), use_container_width=True, hide_index=True)
with tab_match:
    st.dataframe(simple_result_view(matched), use_container_width=True, hide_index=True)
with tab_shift:
    st.dataframe(simple_result_view(possible), use_container_width=True, hide_index=True)
    st.caption("These are putative annotations based on mass-shift agreement, not confirmed PTM/adduct identifications.")
with tab_unassigned:
    st.dataframe(simple_result_view(unassigned), use_container_width=True, hide_index=True)

st.subheader("6. Download")
excel = to_excel_bytes(results, metas, theory)
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.download_button("Download combined Excel", data=excel, file_name="proteasome_mass_analysis.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
with c2:
    st.download_button("Download matched CSV", data=csv_bytes(matched), file_name="matched_masses.csv", mime="text/csv")
with c3:
    st.download_button("Download possible PTM/adduct CSV", data=csv_bytes(possible), file_name="possible_ptm_adduct_masses.csv", mime="text/csv")
with c4:
    st.download_button("Download unassigned CSV", data=csv_bytes(unassigned), file_name="unassigned_masses.csv", mime="text/csv")
