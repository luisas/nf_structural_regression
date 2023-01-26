import pandas as pd
import sys

hits=sys.argv[1]
output=sys.argv[2]
template="template.txt"
output_ids="ids_to_download.txt"
output_chains="chains_summary.txt"


def get_best_hits(hits):

    df = pd.read_csv(hits, sep='\t', header = None)
    # 1 - Get all the sequences with maximum identity
    df_grouped = df.groupby([0]).agg({2:'max'})
    df_grouped = df_grouped.reset_index()
    df_grouped = df_grouped.rename(columns={2:'identity_max'})
    df = pd.merge(df, df_grouped, how='left', on=[0])
    df = df[df[2] == df['identity_max']]

    # Hits presenting a best match that also have the same id name are prioritized
    df["target_id_nochain"] = df[1].str.split("_",expand = True)[0]
    df["target_id_chainmerged"] = df[1].str.replace("_", "").str.lower()

    df_id_match_1 = df[df[0] == df["target_id_nochain"]]
    df_id_match_2 = df[df[0] == df["target_id_chainmerged"]]
    df_id_match = pd.concat([df_id_match_1,df_id_match_2])
    
    # If there are multiple matches per sequence( multiple chains for example), pick the firt one
    df_id_match = df_id_match.groupby(0).first().reset_index()

    # Only retain the dataframe
    df_noid_match = df[~df[0].isin(df_id_match[0])].reset_index(drop = True)
    if(not df_noid_match.empty):
        df_noid_match_filtered = df_noid_match.iloc[[df_noid_match[2].idxmax()]]
        final_df = pd.concat([df_id_match,df_noid_match_filtered])
    else:
        final_df = df_id_match

    return(final_df)

def main():
    df = get_best_hits(hits)
    df.to_csv(output, sep="\t", header=None, index=False)

    # 2. Create file with IDs to download
    df["target_id_nochain"].to_csv(output_ids, sep="\t", header=None, index=False)

    # 3. Create file also with chain informations
    df["chain"] = df[1].str.split("_", expand = True)[1]
    df[[0,1,"chain","target_id_nochain"]].to_csv(output_chains, sep="\t", header=None, index=False)

    # 4. Create template file
    df["sep"] = "_P_"
    df["query_id"]=">"+df[0]
    cols =  ["query_id","sep","target_id_nochain"]
    df[cols].to_csv(template, sep=" ", header=None, index=False)

if __name__ == "__main__":
    main()
