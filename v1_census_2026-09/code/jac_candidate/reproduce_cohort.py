from pathlib import Path
import pandas as pd,json
p=Path(__file__).resolve().parents[2]
out=p/'data/results/jac_candidate';out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(p/'data/raw/emdb_detail.csv',dtype={'deposition':str});c=pd.read_csv(p/'data/raw/emdb_census.csv',dtype={'deposition':str});f=d.resolution_method.fillna('').str.contains('0.143');frame=c[c.deposition.notna()&c.resolution.notna()];post=d.deposition.ge('20220225');oldpost=d.deposition.ge('20220201');early=d.deposition.between('20220201','20220224')
res={'sampling_boundary':'2022-02-01','policy_effective':'2022-02-25','policy_source':'https://www.rcsb.org/news/6218da3152988f064bf8c4a3','sample_n':len(d),'later_stratum_n':int(oldpost.sum()),'later_stratum_non_fsc0143':int((oldpost&~f).sum()),'on_or_after_policy_n':int(post.sum()),'on_or_after_policy_non_fsc0143':int((post&~f).sum()),'between_boundaries_n':int(early.sum()),'between_boundaries_non_fsc0143':d[early&~f].emdb_id.tolist(),'sampling_frame_n':len(frame),'frame_before_sampling_boundary':int(frame.deposition.lt('20220201').sum()),'frame_after_sampling_boundary':int(frame.deposition.ge('20220201').sum())}
w=res['frame_before_sampling_boundary']/len(frame);res['frame_weighted_fsc0143']=float(w*f[~oldpost].mean()+(1-w)*f[oldpost].mean());res['date_subset_methods']={str(k):{'n':int(len(g)),'non_fsc0143':int((~g.resolution_method.fillna('').str.contains('0.143')).sum())} for k,g in d[post].groupby('method')}
assert (res['on_or_after_policy_n'],res['on_or_after_policy_non_fsc0143'])==(781,33);assert res['between_boundaries_non_fsc0143']==['EMD-14430'];(out/'cohort_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
