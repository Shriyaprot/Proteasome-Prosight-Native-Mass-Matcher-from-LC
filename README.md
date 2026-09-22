# Proteasome ProSight Mass Matcher

A Streamlit app for batch analysis of **ProSight Native Interactive Deconvolution DOCX reports** from LC–MS experiments on proteasome subunits.

## What the app reads automatically

For every ProSight Native `.docx` report, the app extracts:

- Word report filename
- LC peak number from names such as `Peak1`, `Peak2`, `Peak3`, ...
- age/sample group (Old/Adult or Infant/Young) when present in the report or RAW filename
- organ (Liver, Lung, Brain, Kidney, Heart, Spleen)
- ProSight **Data File Name** / source RAW file
- retention-time window
- scan numbers
- deconvolution method
- detected neutral mass
- charge-state range
- sum intensity
- relative intensity
- ProSight Top ID fields, if present

No intensity threshold is applied in this version; all reported masses are retained.

## Batch input options

The app supports three ways to analyze many reports:

1. **Upload DOCX files** — select many ProSight reports at once.
2. **Upload a ZIP folder** — zip the experiment folder and upload one archive.
3. **Read a local folder** — enter a folder path when Streamlit is running locally on the same computer.

> A Streamlit Community Cloud deployment cannot directly browse folders on your laptop. Use multi-file upload or ZIP upload there.

## Mass matching

The bundled theoretical database is based on the calculated masses in the supplied mouse 20S proteasome mass list.

For each detected mass, the app finds the nearest theoretical mass and calculates:

- mass difference in Da
- mass error in ppm

A result is assigned only if it is within the user-selected tolerance. Otherwise it remains:

`Unmatched / possible new mass`

The nearest theoretical candidate is still recorded for context, but the app does **not** force-assign an unmatched mass.

## Per-report display

Each DOCX report appears in its own expandable section with:

- metadata and source RAW filename
- a result table
- a graph of **detected mass vs sum intensity**
- matched/unmatched status in the graph hover information
- charge-state and relative-intensity information

This lets each LC elution peak remain traceable to its original report.

## Downloads

The app provides:

- complete Excel workbook
- all masses CSV
- matched/identified masses CSV
- unmatched/possible new masses CSV

The Excel workbook contains separate sheets for:

- `All masses`
- `Matched_identified`
- `Unmatched_new_masses`
- `Report summary`
- `Theory used`

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy from GitHub

1. Create a GitHub repository.
2. Upload the contents of this folder.
3. Commit and push.
4. In Streamlit Community Cloud, create an app from the repository and choose `app.py` as the entry point.

No database or secrets are required.
