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
    "Subunit",
    "Accession",
    "Modification",
    "Found_in_reference",
    "Calculated_Mass_Da",
    "Reference_Measured_Mass_Da",
    "Detected_Mass_Da",
    "Error_STD_Da",
    "Difference_Da",
    "Difference_in_STD",
    "Category",
    "Possible_Explanation",
    "Charge_State_Range",
    "Sum_Intensity",
    "Relative_Intensity_pct",
]


def read_theory(uploaded_file=None) -> pd.DataFrame:
    """Read the reference table while preserving the empirical measured mass and STD.

    The user's Mouse 20S table contains both calculated and historically measured intact
    masses.  Direct matching uses the empirical Measured Mass (Da) and Error (STD) when
    available; Calculated Mass (Da) is retained for reporting.
    """
    if uploaded_file is None:
        theory = pd.read_csv(DEFAULT_THEORY)
    else:
        name = uploaded_file.name.lower()
        theory = pd.read_csv(uploaded_file) if name.endswith(".csv") else pd.read_excel(uploaded_file)

    aliases = {
        "Calculated Mass (Da)": "Calculated_Mass_Da",
        "Calculated_Mass_Da": "Calculated_Mass_Da",
        "Theoretical Mass (Da)": "Calculated_Mass_Da",
        "Theoretical_Mass_Da": "Calculated_Mass_Da",
        "Measured Mass (Da)": "Reference_Measured_Mass_Da",
        "Measured_Mass_Da": "Reference_Measured_Mass_Da",
        "Reference_Measured_Mass_Da": "Reference_Measured_Mass_Da",
        "Error (STD)": "Error_STD_Da",
        "Error_STD": "Error_STD_Da",
        "Error_STD_Da": "Error_STD_Da",
    }
    theory = theory.rename(columns=aliases)

    # Handle CSVs that may contain two differently named "Found in" columns.
    found_cols = [c for c in theory.columns if str(c).lower().startswith("found")]
    if "Found_in_reference" not in theory.columns:
        if found_cols:
            parts = []
            for c in found_cols:
                parts.append(theory[c].fillna("").astype(str).str.strip())
            found = parts[0]
            for p in parts[1:]:
                found = found.where(p.eq(""), found.where(found.eq(""), found + "; ") + p)
            theory["Found_in_reference"] = found
        else:
            theory["Found_in_reference"] = ""

    required = ["Subunit", "Accession", "Modification", "Calculated_Mass_Da"]
    missing = [c for c in required if c not in theory.columns]
    if missing:
        raise ValueError("The theoretical list is missing: " + ", ".join(missing))

    for c in ["Subunit", "Accession"]:
        theory[c] = theory[c].ffill().fillna("")
    theory["Modification"] = theory["Modification"].fillna("")
    theory["Found_in_reference"] = theory["Found_in_reference"].fillna("")

    for c in ["Calculated_Mass_Da", "Reference_Measured_Mass_Da", "Error_STD_Da"]:
        if c not in theory.columns:
            theory[c] = np.nan
        theory[c] = pd.to_numeric(theory[c], errors="coerce")

    keep = [
        "Subunit", "Accession", "Modification", "Found_in_reference",
        "Calculated_Mass_Da", "Reference_Measured_Mass_Da", "Error_STD_Da",
    ]
    return theory[keep].dropna(subset=["Calculated_Mass_Da"]).reset_index(drop=True)


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


def match_masses(df: pd.DataFrame, theory: pd.DataFrame, n_std: float, shift_ppm: float) -> pd.DataFrame:
    """Match intact masses to the empirical reference mass table.

    Primary rule: when a reference measured mass and STD are available, a detected mass
    is a direct match when |detected - measured| <= n_std * STD.  This mirrors the
    empirical tolerance already present in the user's reference workbook.

    If empirical measured/STD values are unavailable for a row, the calculated mass is
    used only as a nearest reference, not as an automatic direct assignment.

    Only masses that fail the direct rule are tested for common PTM/adduct shifts.
    """
    if df.empty:
        return df.copy()

    t = theory.reset_index(drop=True)
    shifts = COMMON_SHIFTS.reset_index(drop=True)
    rows = []

    for _, row in df.iterrows():
        detected = float(row["Detected_Mass_Da"])
        rec = row.to_dict()

        # Direct empirical match: choose the smallest distance in STD units.
        candidates = []
        for ti, ref in t.iterrows():
            measured = ref.get("Reference_Measured_Mass_Da", np.nan)
            std = ref.get("Error_STD_Da", np.nan)
            if pd.notna(measured) and pd.notna(std) and float(std) > 0:
                diff = detected - float(measured)
                z = abs(diff) / float(std)
                candidates.append((z, abs(diff), ti, diff))

        candidates.sort(key=lambda x: (x[0], x[1]))
        direct = candidates[0] if candidates else None

        # Nearest calculated mass is retained for possible secondary interpretation.
        calc = t["Calculated_Mass_Da"].to_numpy(dtype=float)
        nearest_idx = int(np.argmin(np.abs(detected - calc)))
        nearest = t.iloc[nearest_idx]

        rec.update({
            "Subunit": nearest.get("Subunit", ""),
            "Accession": nearest.get("Accession", ""),
            "Modification": nearest.get("Modification", ""),
            "Found_in_reference": nearest.get("Found_in_reference", ""),
            "Calculated_Mass_Da": float(nearest.get("Calculated_Mass_Da")),
            "Reference_Measured_Mass_Da": nearest.get("Reference_Measured_Mass_Da", np.nan),
            "Error_STD_Da": nearest.get("Error_STD_Da", np.nan),
            "Difference_Da": detected - float(nearest.get("Reference_Measured_Mass_Da")) if pd.notna(nearest.get("Reference_Measured_Mass_Da", np.nan)) else detected - float(nearest.get("Calculated_Mass_Da")),
            "Difference_in_STD": np.nan,
            "Category": "Unassigned / investigate",
            "Possible_Explanation": "",
        })

        if direct is not None and direct[0] <= n_std:
            z_abs, _, ti, diff = direct
            ref = t.iloc[ti]
            std = float(ref["Error_STD_Da"])
            rec.update({
                "Subunit": ref.get("Subunit", ""),
                "Accession": ref.get("Accession", ""),
                "Modification": ref.get("Modification", ""),
                "Found_in_reference": ref.get("Found_in_reference", ""),
                "Calculated_Mass_Da": float(ref["Calculated_Mass_Da"]),
                "Reference_Measured_Mass_Da": float(ref["Reference_Measured_Mass_Da"]),
                "Error_STD_Da": std,
                "Difference_Da": diff,
                "Difference_in_STD": diff / std,
                "Category": "Matched",
                "Possible_Explanation": "",
            })
            rows.append(rec)
            continue

        # Secondary PTM/adduct annotation: compare shifts from each known reference mass.
        best = None
        for ti, ref in t.iterrows():
            base = ref.get("Reference_Measured_Mass_Da", np.nan)
            if pd.isna(base):
                base = ref.get("Calculated_Mass_Da", np.nan)
            if pd.isna(base):
                continue
            base = float(base)
            for _, s in shifts.iterrows():
                shift_da = float(s["Shift_Da"])
                for direction in (1.0, -1.0):
                    expected = base + direction * shift_da
                    resid_da = detected - expected
                    resid_ppm = ppm_error(detected, expected)
                    score = abs(resid_ppm)
                    if best is None or score < best[0]:
                        best = (score, ti, s, direction, resid_da, resid_ppm, expected)

        if best is not None and best[0] <= shift_ppm:
            _, ti, s, direction, resid_da, resid_ppm, expected = best
            ref = t.iloc[ti]
            sign = "+" if direction > 0 else "−"
            category = "Possible PTM" if s["Shift_Type"] == "PTM" else "Possible adduct"
            base = ref.get("Reference_Measured_Mass_Da", np.nan)
            if pd.isna(base):
                base = ref.get("Calculated_Mass_Da", np.nan)
            rec.update({
                "Subunit": ref.get("Subunit", ""),
                "Accession": ref.get("Accession", ""),
                "Modification": ref.get("Modification", ""),
                "Found_in_reference": ref.get("Found_in_reference", ""),
                "Calculated_Mass_Da": float(ref["Calculated_Mass_Da"]),
                "Reference_Measured_Mass_Da": ref.get("Reference_Measured_Mass_Da", np.nan),
                "Error_STD_Da": ref.get("Error_STD_Da", np.nan),
                "Difference_Da": detected - float(base),
                "Difference_in_STD": (detected - float(base)) / float(ref["Error_STD_Da"]) if pd.notna(ref.get("Error_STD_Da", np.nan)) and float(ref["Error_STD_Da"]) > 0 else np.nan,
                "Category": category,
                "Possible_Explanation": f"{sign}{s['Shift_Name']} ({sign}{shift_da:.4f} Da); residual {resid_da:+.4f} Da ({resid_ppm:+.2f} ppm)",
            })

        rows.append(rec)

    return pd.DataFrame(rows)


def simple_result_view(result: pd.DataFrame) -> pd.DataFrame:
    return result[[c for c in RESULT_COLUMNS if c in result.columns]].copy()


def report_plot(result: pd.DataFrame):
    p = result.copy()
    p["Assignment"] = np.where(
        p["Category"].eq("Matched"),
        p["Subunit"].fillna("") + " — " + p["Modification"].fillna(""),
        np.where(
            p["Category"].isin(["Possible PTM", "Possible adduct"]),
            p["Category"] + ": " + p["Possible_Explanation"].fillna(""),
            "Unassigned",
        ),
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
    """Create ONE simple combined worksheet, grouped by report.

    Report metadata is written once above each report's mass table.  No extra sheets are
    created for categories, theory, or shift references.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Combined results"
    ws.sheet_view.showGridLines = False

    title_fill = PatternFill("solid", fgColor="1F4E78")
    title_font = Font(color="FFFFFF", bold=True, size=12)
    meta_fill = PatternFill("solid", fgColor="EAF2F8")

    row = 1
    for result, meta in zip(results, metas):
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=max(8, len(RESULT_COLUMNS)))
        tc = ws.cell(row, 1, f"{meta.get('Report_File','')}  |  {meta.get('Peak','')}")
        tc.fill = title_fill
        tc.font = title_font
        row += 1

        # Keep only the essential source identifiers at the top of each report block.
        for label, value in [
            ("Report file", meta.get("Report_File", "")),
            ("Data file", meta.get("Data_File_Name", "")),
        ]:
            ws.cell(row, 1, label).font = Font(bold=True)
            ws.cell(row, 1).fill = meta_fill
            ws.cell(row, 2, value)
            row += 1

        row += 1
        row = append_result_table(ws, row, result)
        row += 2

    for idx in range(1, ws.max_column + 1):
        max_len = 0
        for cell in ws[get_column_letter(idx)]:
            if cell.value is not None:
                max_len = max(max_len, min(len(str(cell.value)), 48))
        ws.column_dimensions[get_column_letter(idx)].width = max(12, min(max_len + 2, 48))

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
    n_std = st.number_input("Direct match window (× reference STD)", min_value=0.5, value=2.0, step=0.5)
    shift_ppm = st.number_input("PTM/adduct residual tolerance (ppm)", min_value=0.1, value=5.0, step=0.5)
    st.caption("Direct matching uses the reference Measured Mass ± N×Error(STD). Default: ±2 STD. PTM/adduct checks are only attempted after direct matching fails.")

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
        result = match_masses(parsed, theory, n_std, shift_ppm)
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
matched = combined[combined["Category"] == "Matched"].copy()
possible = combined[combined["Category"].isin(["Possible PTM", "Possible adduct"])].copy()
unassigned = combined[combined["Category"] == "Unassigned / investigate"].copy()
st.dataframe(simple_result_view(combined), use_container_width=True, hide_index=True)
st.caption("Possible PTM/adduct labels are putative mass-shift annotations, not confirmed identities.")

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
