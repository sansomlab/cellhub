import pandas as pd
import numpy as np

def preprocess_cellranger_stats(infile, outfile):
    libraries = pd.read_csv(infile, sep='\t')
    libraries.set_index("library_id", inplace=True)

    ss = []

    for library_id in libraries.index:
        cellranger_summary = "/".join(["cellranger.multi.dir", library_id,
                                       "outs/per_sample_outs", library_id,
                                       "metrics_summary.csv"])
        map_sum = pd.read_csv(cellranger_summary, sep=',')

        # (need the group name to ensure entries are unique)
        map_sum = map_sum.replace(np.nan, 'na', regex=True)
        map_sum["id"] =  map_sum["Library or Sample"] + \
                         "_" + map_sum["Library Type"] + \
                         "_" + map_sum["Group Name"] + \
                         "_" + map_sum["Metric Name"]

        map_sum["id"] = [ x.replace(' ','_').lower() for x in map_sum["id"].values]
        map_sum.index = map_sum["id"]

        sub_map_sum = pd.DataFrame(map_sum["Metric Value"]).transpose()

        col_names = []
        for col in sub_map_sum.columns:

            x = sub_map_sum[col].values[0]

            if x[-1:]=="%":
                col_names.append(col + "_pct")
                x = x.replace("%","")
            else:
                col_names.append(col)

            x = x.replace(",","")


            sub_map_sum[col] = pd.to_numeric(x)
        sub_map_sum.columns = col_names
        sub_map_sum['library_id'] = library_id

        ss.append(sub_map_sum)

    summary = pd.concat(ss)
    summary.to_csv(outfile, sep = '\t', index=False)
