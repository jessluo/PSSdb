# Pelagic Size Structure database (PSSdb)
A collaborative project between NOAA/Geophysical Fluid Dynamics Laboratory and Sorbonne Université/Laboratoire d’Océanographie de Villefranche-sur-Mer (LOV)

v1.2: step4 has been modified to calculate PSD and save intermediate data products of NBSS and PSD calculations

### Repository organization

```
PSSdb/
  .gitignore
  requirements.txt
  raw/
    project_list_all.xslx: list of projects hosted on Ecotaxa, Ecopart, and IFCB dashboards
    project_IFCB/UVP/Zooscan_standardizer.xlsx: spreadsheets to map and standardize exported projects
    ecotaxa/
           instrument/ecotaxa_export_projectID_exportdateGMT_exporttimeGMT.tsv: native export file for projects hosted on Ecotaxa
    ecopart/
           instrument/ecopart_export_raw/detailed_projectID_exportdateGMT_exporttimeGMT.tsv: native export file for projects hosted on Ecopart
    flags/
           instrument/project_projectID_flags.tsv: table including flagged samples based on 5 criteria and overruling test
    IFCB_dashboard_downloads/
    raw_standardized/
           instrument/standardized_project_projectID.tsv: standardized table
  reports/
    instrument/report_project_projectID.html: an interactive report showing the flagged samples locations, number of ROI, percentage of validation, and flag details
  scripts/
    step0_list_projects.py
    step1_export_projects.py
    step2_standardize_projects.py
    step3_grid_projects.py
    step4_compute_NBSS.py
```
        
### Hints on code documentation

https://realpython.com/documenting-python-code/
