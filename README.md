# Proteasome ProSight Native Mass Matcher

Streamlit app for batch parsing ProSight Native deconvolution DOCX reports and matching intact proteasome masses to a reference table.

## V4 matching logic

- Direct match uses the empirical `Measured Mass (Da)` and `Error (STD)` from the reference table.
- Default direct-match window is +/- 2 x STD and is adjustable in the Streamlit sidebar.
- `Calculated Mass (Da)` is retained next to the empirical reference and detected mass for interpretation.
- PTM/adduct checks are attempted only after direct matching fails.
- Common PTM/adduct suggestions are putative annotations, not confirmed identifications.

## Combined Excel output

The downloaded workbook contains one sheet only: `Combined results`.
For each ProSight report, `Report file` and `Data file` appear once, followed by a concise mass table. Charge state and intensity columns are placed after the main identification/mass-matching columns.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```
