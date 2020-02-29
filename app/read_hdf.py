import pandas as pd

def read_cag_details(hdf_fp, cag_id):
  return pd.read_hdf(
    hdf_fp,
    "/annot/gene/all",
    where = "CAG == {}".format(cag_id)
  ).fillna("").applymap(str)

def read_cag_abundance(hdf_fp, cag_id):
  return pd.DataFrame({
    "abundance": pd.read_hdf(
      hdf_fp,
      "/abund/cag/wide",
      where = "CAG == {}".format(cag_id)
    ).set_index(
      "CAG"
    ).T[cag_id]
  }).reset_index(
  ).rename(
    columns = {
      "index": "specimen"
    }
  )
  

def read_multiple_cag_abundances(hdf_fp, cag_id_list):
  if isinstance(cag_id_list, list) is False:
    cag_id_list = [cag_id_list]

  return pd.DataFrame({
    str(cag_id): pd.read_hdf(
      hdf_fp,
      "/abund/cag/wide",
      where = "CAG == {}".format(cag_id)
    ).set_index(
      "CAG"
    ).T[int(cag_id)]
    for cag_id in cag_id_list
  }).reset_index(
  ).rename(
    columns = {
      "index": "specimen"
    }
  )
