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

# Putative intact-mass shifts for SECONDARY annotation only.
# A shift match is a hypothesis, not a confirmed PTM/adduct identity.
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
    "Calculated_Mass_Da",
    "Detected_Mass_Da",
    "Difference_Da",
    "Difference_ppm",
    "Category",
    "Charge_State_Range",
    "Charge_State_Count",
    "Sum_Intensity",
    "Relative_Intensity_pct",
    "Possible_Explanation",
]



MATCHED_COLUMNS = [
    "Subunit",
    "Accession",
    "Modification",
    "Calculated_Mass_Da",
    "Occurrence_Count",
    "Detected_Mass_Da",
    "Difference_Da",
    "Difference_ppm",
    "Charge_State_Range",
    "Charge_State_Count",
    "Sum_Intensity",
    "Relative_Intensity_pct",
    "Age_Group",
    "Organ",
    "Peak",
    "Report_File",
    "Data_File_Name",
]

TRACKER_COLUMNS = [
    "Cluster_Mass_Da",
    "Occurrence_Count",
    "Detected_Mass_Da",
    "Age_Group",
    "Organ",
    "Peak",
    "Report_File",
    "Data_File_Name",
    "Category",
    "Subunit",
    "Accession",
    "Modification",
    "Calculated_Mass_Da",
    "Difference_Da",
    "Difference_ppm",
    "Charge_State_Range",
    "Charge_State_Count",
    "Sum_Intensity",
    "Relative_Intensity_pct",
    "Possible_Explanation",
]



def read_theory(uploaded_file=None) -> pd.DataFrame:
    """Read the theoretical proteasome mass list.

    Only calculated/theoretical masses are used. Historical measured masses,
    experiment-specific STD values, and previous-study tissue annotations are ignored.
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
    }
    theory = theory.rename(columns=aliases)

    required = ["Subunit", "Accession", "Modification", "Calculated_Mass_Da"]
    missing = [c for c in required if c not in theory.columns]
    if missing:
        raise ValueError("The theoretical list is missing: " + ", ".join(missing))

    for c in ["Subunit", "Accession"]:
        theory[c] = theory[c].ffill().fillna("")
    theory["Modification"] = theory["Modification"].fillna("")
    theory["Calculated_Mass_Da"] = pd.to_numeric(theory["Calculated_Mass_Da"], errors="coerce")

    keep = ["Subunit", "Accession", "Modification", "Calculated_Mass_Da"]
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
            "Charge_State_Count": (hi - lo + 1) if pd.notna(lo) and pd.notna(hi) else np.nan,
            "Sum_Intensity": pd.to_numeric(block.get("Sum Intensity", ""), errors="coerce"),
            "Relative_Intensity_pct": pd.to_numeric(
                str(block.get("Relative Intensity", "")).replace("%", ""), errors="coerce"
            ),
        })
        i = max(j, i + 1)

    return pd.DataFrame(records), meta


def ppm_error(observed: float, expected: float) -> float:
    return (observed - expected) / expected * 1_000_000.0


def match_masses(
    df: pd.DataFrame,
    theory: pd.DataFrame,
    tolerance_mode: str,
    direct_tolerance: float,
    near_tolerance: float,
    shift_ppm: float,
) -> pd.DataFrame:
    """Match detected masses against calculated theoretical masses only.

    No historical experimental STD is used. Direct matching is controlled by the user
    in either Da or ppm. A second, wider near-theoretical window can be used to retain
    borderline candidates for recurrence tracking. PTM/adduct annotation is secondary.
    """
    if df.empty:
        return df.copy()

    t = theory.reset_index(drop=True)
    shifts = COMMON_SHIFTS.reset_index(drop=True)
    rows = []

    def metric(abs_da: float, calc: float) -> float:
        if tolerance_mode == "ppm":
            return abs_da / calc * 1_000_000.0
        return abs_da

    for _, row in df.iterrows():
        detected = float(row["Detected_Mass_Da"])
        rec = row.to_dict()

        candidates = []
        for ti, ref in t.iterrows():
            calc = float(ref["Calculated_Mass_Da"])
            diff = detected - calc
            score = metric(abs(diff), calc)
            candidates.append((score, abs(diff), ti, diff))

        candidates.sort(key=lambda x: (x[0], x[1]))
        score, _, nearest_ti, diff = candidates[0]
        nearest = t.iloc[nearest_ti]
        calc = float(nearest["Calculated_Mass_Da"])
        diff_ppm = ppm_error(detected, calc)

        rec.update({
            "Subunit": nearest.get("Subunit", ""),
            "Accession": nearest.get("Accession", ""),
            "Modification": nearest.get("Modification", ""),
            "Calculated_Mass_Da": calc,
            "Difference_Da": diff,
            "Difference_ppm": diff_ppm,
            "Category": "Unassigned / investigate",
            "Possible_Explanation": "",
        })

        if score <= direct_tolerance:
            rec["Category"] = "Matched"
            rows.append(rec)
            continue

        if score <= near_tolerance:
            rec["Category"] = "Near theoretical / investigate"

        # Secondary PTM/adduct hypothesis against each theoretical proteoform.
        best = None
        for ti, ref in t.iterrows():
            base = float(ref["Calculated_Mass_Da"])
            for _, s in shifts.iterrows():
                shift_da = float(s["Shift_Da"])
                for direction in (1.0, -1.0):
                    expected = base + direction * shift_da
                    resid_da = detected - expected
                    resid_ppm = ppm_error(detected, expected)
                    score_shift = abs(resid_ppm)
                    if best is None or score_shift < best[0]:
                        best = (score_shift, ti, s, direction, resid_da, resid_ppm)

        if best is not None and best[0] <= shift_ppm:
            _, ti, s, direction, resid_da, resid_ppm = best
            ref = t.iloc[ti]
            base = float(ref["Calculated_Mass_Da"])
            sign = "+" if direction > 0 else "−"
            category = "Possible PTM" if s["Shift_Type"] == "PTM" else "Possible adduct"
            rec.update({
                "Subunit": ref.get("Subunit", ""),
                "Accession": ref.get("Accession", ""),
                "Modification": ref.get("Modification", ""),
                "Calculated_Mass_Da": base,
                "Difference_Da": detected - base,
                "Difference_ppm": ppm_error(detected, base),
                "Category": category,
                "Possible_Explanation": (
                    f"{sign}{s['Shift_Name']} ({sign}{float(s['Shift_Da']):.4f} Da); "
                    f"residual {resid_da:+.4f} Da ({resid_ppm:+.2f} ppm)"
                ),
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
            p["Category"],
        ),
    )
    fig = px.scatter(
        p,
        x="Detected_Mass_Da",
        y="Sum_Intensity",
        symbol="Category",
        size="Relative_Intensity_pct",
        hover_data={
            "Detected_Mass_Da": ":.2f",
            "Sum_Intensity": ":.3g",
            "Relative_Intensity_pct": True,
            "Charge_State_Range": True,
            "Assignment": True,
            "Category": True,
        },
        title="Deconvoluted masses",
        labels={"Detected_Mass_Da": "Detected mass (Da)", "Sum_Intensity": "Sum intensity"},
    )
    fig.update_traces(marker={"line": {"width": 1}})
    fig.update_layout(height=420, legend_title_text="Classification")
    return fig


def add_metadata_columns(result: pd.DataFrame, meta: dict[str, str]) -> pd.DataFrame:
    out = result.copy()
    for col, key in [
        ("Age_Group", "Age_Group"),
        ("Organ", "Organ"),
        ("Peak", "Peak"),
        ("Report_File", "Report_File"),
        ("Data_File_Name", "Data_File_Name"),
    ]:
        out[col] = meta.get(key, "")
    return out



def compile_matched_subunits(results: list[pd.DataFrame], metas: list[dict[str, str]]) -> pd.DataFrame:
    """Compile every direct theoretical match across all uploaded reports.

    One row is kept per matched observation so the user can see where each subunit/proteoform
    was detected and inspect charge-state and intensity evidence. Occurrence_Count gives the
    total number of times that exact Subunit + Accession + Modification was matched in the batch.
    """
    pieces = []
    for result, meta in zip(results, metas):
        matched = result[result["Category"] == "Matched"].copy()
        if matched.empty:
            continue
        matched = add_metadata_columns(matched, meta)
        pieces.append(matched)

    if not pieces:
        return pd.DataFrame(columns=MATCHED_COLUMNS)

    out = pd.concat(pieces, ignore_index=True)
    group_cols = ["Subunit", "Accession", "Modification"]
    counts = out.groupby(group_cols, dropna=False)["Detected_Mass_Da"].transform("size")
    out["Occurrence_Count"] = counts.astype(int)

    keep = [c for c in MATCHED_COLUMNS if c in out.columns]
    out = out[keep].sort_values(
        ["Subunit", "Modification", "Organ", "Age_Group", "Peak", "Detected_Mass_Da"],
        na_position="last",
    ).reset_index(drop=True)
    return out

def cluster_candidate_masses(candidates: pd.DataFrame, tolerance_da: float) -> pd.DataFrame:
    """Group recurring candidate masses by observed mass similarity.

    Clustering is based on detected mass in Da. This keeps recurring candidate masses
    together without assuming an experimental STD that has not been measured in this dataset.
    """
    if candidates.empty:
        return candidates.copy()

    c = candidates.sort_values("Detected_Mass_Da").reset_index(drop=True).copy()
    cluster_ids = []
    centers: list[list[float]] = []

    for mass in c["Detected_Mass_Da"].astype(float):
        if not centers:
            centers.append([mass])
            cluster_ids.append(0)
            continue
        current_center = float(np.median(centers[-1]))
        if abs(mass - current_center) <= tolerance_da:
            centers[-1].append(mass)
            cluster_ids.append(len(centers) - 1)
        else:
            centers.append([mass])
            cluster_ids.append(len(centers) - 1)

    c["_cluster_id"] = cluster_ids
    center_map = {i: float(np.median(vals)) for i, vals in enumerate(centers)}
    count_map = c.groupby("_cluster_id").size().to_dict()
    c["Cluster_Mass_Da"] = c["_cluster_id"].map(center_map)
    c["Occurrence_Count"] = c["_cluster_id"].map(count_map)
    c = c.drop(columns=["_cluster_id"])

    keep = [col for col in TRACKER_COLUMNS if col in c.columns]
    return c[keep].sort_values(
        ["Occurrence_Count", "Cluster_Mass_Da", "Detected_Mass_Da"],
        ascending=[False, True, True],
    ).reset_index(drop=True)


def style_result_header(ws, row: int, headers: list[str]):
    for col, h in enumerate(headers, 1):
        c = ws.cell(row, col, h)
        c.font = Font(bold=True)
        c.fill = PatternFill("solid", fgColor="D9EAF7")
        c.alignment = Alignment(wrap_text=True)


def append_result_table(ws, start_row: int, result: pd.DataFrame) -> int:
    view = simple_result_view(result)
    headers = list(view.columns)
    style_result_header(ws, start_row, headers)
    for r_idx, values in enumerate(view.itertuples(index=False, name=None), start_row + 1):
        for c_idx, value in enumerate(values, 1):
            if pd.isna(value):
                value = None
            cell = ws.cell(r_idx, c_idx, value)
            if headers[c_idx - 1] == "Charge_State_Range":
                cell.number_format = "@"
    return start_row + len(view) + 1


def autosize(ws, max_width: int = 48):
    for idx in range(1, ws.max_column + 1):
        max_len = 0
        for cell in ws[get_column_letter(idx)]:
            if cell.value is not None:
                max_len = max(max_len, min(len(str(cell.value)), max_width))
        ws.column_dimensions[get_column_letter(idx)].width = max(12, min(max_len + 2, max_width))


def write_dataframe_sheet(ws, df: pd.DataFrame):
    headers = list(df.columns)
    style_result_header(ws, 1, headers)
    for r_idx, values in enumerate(df.itertuples(index=False, name=None), 2):
        for c_idx, value in enumerate(values, 1):
            if pd.isna(value):
                value = None
            cell = ws.cell(r_idx, c_idx, value)
            if headers[c_idx - 1] == "Charge_State_Range":
                cell.number_format = "@"
    ws.freeze_panes = "A2"
    ws.sheet_view.showGridLines = False
    autosize(ws)


def to_excel_bytes(
    results: list[pd.DataFrame],
    metas: list[dict[str, str]],
    candidate_tracker: pd.DataFrame,
    matched_subunits: pd.DataFrame,
) -> bytes:
    """Create one workbook with combined results, candidate tracking, and compiled matches."""
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

    autosize(ws)

    tracker_ws = wb.create_sheet("Encountered candidates")
    if candidate_tracker.empty:
        tracker_ws["A1"] = "No unmatched / near-theoretical / possible PTM-adduct masses were found."
    else:
        write_dataframe_sheet(tracker_ws, candidate_tracker)

    matched_ws = wb.create_sheet("Matched subunits")
    if matched_subunits.empty:
        matched_ws["A1"] = "No theoretical proteasome masses were matched in this batch."
    else:
        write_dataframe_sheet(matched_ws, matched_subunits)

    output = io.BytesIO()
    wb.save(output)
    return output.getvalue()


def combined_csv_bytes(df: pd.DataFrame) -> bytes:
    safe = simple_result_view(df)
    return safe.to_csv(index=False).encode("utf-8-sig")


def load_docx_sources(mode: str):
    items: list[tuple[str, bytes]] = []
    if mode == "Upload DOCX files":
        uploaded = st.file_uploader(
            "Select all ProSight Native DOCX reports", type=["docx"], accept_multiple_files=True
        )
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
        folder = st.text_input(
            "Local folder path", placeholder=r"C:\Users\you\Documents\ProSight_reports"
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
                        items.append((path.name, path.read_bytes()))
    return items


st.set_page_config(page_title="Proteasome ProSight Mass Matcher", page_icon="🧬", layout="wide")
st.title("🧬 Proteasome ProSight Mass Matcher")
st.caption(
    "Batch-read ProSight Native DOCX reports, compare detected intact masses with calculated "
    "proteasome masses, retain intensity/charge-state evidence, flag putative PTM/adduct shifts, "
    "and track recurring candidate masses."
)

with st.sidebar:
    st.header("Matching settings")
    tolerance_mode_label = st.radio(
        "Direct matching tolerance", ["Da", "ppm"], horizontal=True,
        help="No historical STD is used. Choose an absolute Da window or a ppm window."
    )
    tolerance_mode = "Da" if tolerance_mode_label == "Da" else "ppm"
    if tolerance_mode == "Da":
        direct_tolerance = st.number_input(
            "Direct match window (± Da)", min_value=0.1, max_value=50.0, value=3.0, step=0.5
        )
        near_tolerance = st.number_input(
            "Track near-theoretical masses up to (± Da)",
            min_value=float(direct_tolerance), max_value=100.0,
            value=max(5.0, float(direct_tolerance)), step=0.5
        )
    else:
        direct_tolerance = st.number_input(
            "Direct match window (± ppm)", min_value=1.0, max_value=1000.0, value=100.0, step=10.0
        )
        near_tolerance = st.number_input(
            "Track near-theoretical masses up to (± ppm)",
            min_value=float(direct_tolerance), max_value=5000.0,
            value=max(200.0, float(direct_tolerance)), step=10.0
        )
    shift_ppm = st.number_input(
        "PTM/adduct residual tolerance (ppm)", min_value=0.1, value=5.0, step=0.5
    )
    recurrence_da = st.number_input(
        "Recurring candidate grouping tolerance (Da)", min_value=0.1, value=3.0, step=0.5
    )
    st.caption(
        "Direct matching uses calculated/theoretical masses only. Historical measured masses and STD values are not used. "
        "Intensity and charge-state evidence are shown next to every assignment so you can judge whether the detected peak is convincing."
    )

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
mode = st.radio(
    "How do you want to provide the reports?",
    ["Upload DOCX files", "Upload a ZIP folder", "Read a local folder"],
    horizontal=True,
)
if mode == "Read a local folder":
    st.info(
        "Local-folder mode works only when Streamlit is running on the same computer as the folder. "
        "Streamlit Community Cloud cannot browse a folder on your laptop."
    )
sources = load_docx_sources(mode)
if not sources:
    st.info("Add ProSight DOCX reports to start the analysis.")
    st.stop()

results: list[pd.DataFrame] = []
metas: list[dict[str, str]] = []
tracker_inputs: list[pd.DataFrame] = []
parse_errors = []

for filename, raw in sources:
    try:
        parsed, meta = parse_prosight_docx(raw, filename)
        if parsed.empty:
            parse_errors.append(f"{filename}: no deconvoluted mass blocks were found")
            continue
        result = match_masses(parsed, theory, tolerance_mode, direct_tolerance, near_tolerance, shift_ppm)
        results.append(result)
        metas.append(meta)

        nonmatched = result[result["Category"] != "Matched"].copy()
        if not nonmatched.empty:
            tracker_inputs.append(add_metadata_columns(nonmatched, meta))
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
tracker_raw = pd.concat(tracker_inputs, ignore_index=True) if tracker_inputs else pd.DataFrame()
candidate_tracker = cluster_candidate_masses(tracker_raw, recurrence_da)
matched_subunits = compile_matched_subunits(results, metas)

st.subheader("3. Batch overview")
counts = combined["Category"].value_counts()
a, b, c, d, e = st.columns(5)
a.metric("DOCX reports", len(results))
b.metric("Detected masses", len(combined))
c.metric("Matched", int(counts.get("Matched", 0)))
d.metric(
    "Possible PTM/adduct",
    int(counts.get("Possible PTM", 0) + counts.get("Possible adduct", 0)),
)
e.metric("Tracked candidates", len(candidate_tracker))

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
st.dataframe(simple_result_view(combined), use_container_width=True, hide_index=True)
st.caption(
    "Matched = within the selected Da/ppm window of a calculated theoretical mass. "
    "Possible PTM/adduct annotations are mass-shift hypotheses only."
)

if not candidate_tracker.empty:
    with st.expander("Encountered candidate masses across files"):
        st.dataframe(candidate_tracker, use_container_width=True, hide_index=True)

if not matched_subunits.empty:
    with st.expander("Compiled matched subunits and modifications"):
        st.dataframe(matched_subunits, use_container_width=True, hide_index=True)

st.subheader("6. Download")
excel = to_excel_bytes(results, metas, candidate_tracker, matched_subunits)
c1, c2 = st.columns(2)
with c1:
    st.download_button(
        "Download combined Excel",
        data=excel,
        file_name="proteasome_mass_analysis.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
with c2:
    st.download_button(
        "Download combined CSV",
        data=combined_csv_bytes(combined),
        file_name="proteasome_mass_analysis.csv",
        mime="text/csv",
    )
