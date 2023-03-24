from glob import glob
import os, json, argparse
from datetime import datetime

def convert_drugname(abbr):

    drug_mapper = {"AMC": "AMOXICILIN-CLAVULANATE","AMI": "AMIKACIN","AMX": "AMOXICILIN","AZM": "AZITHROMYCIN","BDQ": "BEDAQUILINE","CAP": "CAPREOMYCIN",
                   "CFZ": "CLOFAZIMINE","CIP": "CIPROFLOXACIN","CLR": "CLARITHROMYCIN","CYC": "CYCLOSERINE","DCS": "D-CYCLOSERINE","DLM": "DELAMANID",
                   "EMB": "ETHAMBUTOL","ETH": "ETHIONAMIDE","ETP": "ERTAPENEM","FQS": "FLUOROQUINOLONE","GEN": "GENTAMICIN","GFX": "GATIFLOXACIN",
                   "IMI": "IMIPENEM","INH": "ISONIAZID","KAN": "KANAMYCIN","LEV": "LEVOFLOXACIN","LZD": "LINEZOLID","MEF": "MEFLOQUINE","MPM": "MEROPENEM",
                   "MXF": "MOXIFLOXACIN","OFX": "OFLOXACIN","PAN": "PRETOMANID","PAS": "PAS","PTO": "PROTHIONAMIDE","PZA": "PYRAZINAMIDE","RFB": "RIFABUTIN",
                   "RIF": "RIFAMPICIN","STM": "STREPTOMYCIN","STX": "SITAFLOXACIN","SXT": "COTRIMOXAZOLE","SZD": "SUTEZOLID","TRD": "TERIZIDONE","TZE": "THIOACETAZONE"}
    try:
        return drug_mapper[abbr].lower()
    except KeyError:
        return abbr    

def collect_tbprofilers(rep_dir,report_dict):
    
    # Extract tbprofiler dictionary
    tbp_out_lst = sorted(glob(os.path.join(rep_dir,"*tbprofiler.results.json")))
    catalog= "TB-Profiler"
    for fj in tbp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        if sid not in report_dict:
            report_dict[sid] = {"seqid":sid,"lineage":{},"res_var":{},"other_var":{}, "valid":1}
        ptr = report_dict[sid]
        with open(fj) as hdl:
            res_dict = json.load(hdl)
            # Directly extracting the sublineage
            if "main_lin" in res_dict:
                   pp = res_dict["main_lin"].split(';')
                   if len(pp)==1:
                    if len(res_dict["lineage"]):
                      ptr["lineage"] = res_dict["lineage"][-1] # selecting the last sublineage
                   elif len(pp)>1:
                        ptr["lineage"] = {"lin":res_dict["main_lin"],"family": "Metagenomic sample: multiple lineages were detected.",
                                          "spoligotype": "","rd": "","frac":'' }
                        ptr["valid"]=0

            ptr["res_var"][catalog] = {}
            for rvr in res_dict["dr_variants"]:
                ptr_cat = ptr["res_var"][catalog]
                for drg in rvr["drugs"]:
                    if drg["drug"] not in ptr_cat:
                            ptr_cat[drg["drug"]]= []
                    ptr_cat[drg["drug"]].append({"mutation":f"{rvr['gene']}@{rvr['change']}", "phenotype":"R","freq":round(rvr['freq'],2)})

            if catalog not in ptr["other_var"]:
                ptr["other_var"][catalog] = {}
            ptr_cat = ptr["other_var"][catalog]
            
            for ovr in res_dict["other_variants"]:
                for drg in ovr["gene_associated_drugs"]:
                    if drg not in ptr_cat:
                            ptr_cat[drg]= []
                    ptr_cat[drg].append({"mutation":f"{ovr['gene']}@{ovr['change']}","phenotype":'U',"freq":round(rvr['freq'],2)})

    return report_dict

def collect_cryptics(rep_dir,report_dict):
    cryp_out_lst = glob(os.path.join(rep_dir,"*cryptic.json"))
    
    for fj in cryp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        catalog = os.path.basename(fj).split('.')[1]

        if sid not in report_dict:
            report_dict[sid] = {"seqid":sid,"lineage":{},"res_var":{},"other_var":{}}
        ptr = report_dict[sid]
        with open(fj) as hdl:
            res_dict = json.load(hdl)
            res_dict = res_dict["data"]

        if catalog not in ptr["res_var"]:
            ptr["res_var"][catalog] = {}
        ptr_cat = ptr["res_var"][catalog]

        if "EFFECTS" in res_dict:
            for drg, mut_list in res_dict["EFFECTS"].items():
                drg = convert_drugname(drg)

                for mut in mut_list:
                    
                    if "PHENOTYPE" in mut:
                        continue

                    if mut['PREDICTION'] == "S":
                        if catalog not in ptr["other_var"]:
                            ptr["other_var"][catalog] = {}
                        tptr = ptr["other_var"][catalog]
                        if drg not in tptr:
                            tptr[drg]=[]
                        tptr[drg].append({"mutation":f"{mut['GENE']}@{mut['MUTATION']}","phenotype":'S'})
                    else:
                        if drg not in ptr_cat:
                            ptr_cat[drg]= []
                        ptr_cat[drg].append({"mutation":f"{mut['GENE']}@{mut['MUTATION']}", "phenotype":mut['PREDICTION']})
    
    return report_dict

def generate(in_dir,template):
    report_dict = {}
    collect_tbprofilers(in_dir,report_dict)
    collect_cryptics(in_dir,report_dict)
    djson = json.dumps(report_dict)
    report_lines = []
    with open(template) as hdl:
        for ll in hdl:
            if ll.find('<!--DATASLOT-->')>-1:
                ll = f"<script>let djson={djson}</script>\n";
            report_lines.append(ll)
    dt_tag = datetime.now().strftime("%y%m%d-%H%M%S")
    with open(f"spear-mtb_report_{dt_tag}.html",'w') as hdl:
        hdl.write(''.join(report_lines))
    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir",required=True,help="Iput directory containing results of ")
    parser.add_argument("--html",required=True,help="The report html template")

    options = parser.parse_args()

    generate(options.in_dir,options.html)
    
    
