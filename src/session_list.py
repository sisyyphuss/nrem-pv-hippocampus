# Full chronologically ordered session paths for the 3 subjects retained
# for the primary analysis after the structural QC pass (see notebook
# section 2, "Data QC and session selection"):
#   - DT2: identical task ('wheel_alternation_task') on all 11 days -> cleanest subject
#   - DT5: mostly 'central_alternation_task', a few off-task/zero-unit days (excluded per-pair)
#   - DT7: mostly 'wheel_alternation_task', a few off-task/zero-unit days (excluded per-pair)
# Excluded entirely: DT8 (sleep_states missing on all 17 sessions),
#                     DT9 (task rotates across 4 different apparatuses almost daily,
#                          and last 5/14 sessions have 0 hippocampal units).
SESSIONS = {
"DT2": [
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3324um-1000um-20160217-160217-103121_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3324um-1144um-20160219-160219-130435_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3396um-1180um-20160220-160220-141839_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3396um-1180um-20160221-160221-135947_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3468um-1216um-20160222-160222-120902_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3468um-1216um-20160223-160223-121321_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3468um-1216um-20160225-160225-130701_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3468um-1252um-20160226-160226-124132_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3540um-1288um-20160227-160227-121226_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3540um-1288um-20160228-160228-135001_behavior+ecephys.nwb",
 "sub-DT2/sub-DT2_ses-DT2-rPPC-rCCG-3612um-1360um-20160304-160304-150220_behavior+ecephys.nwb",
],
"DT5": [
 "sub-DT5/sub-DT5_ses-20161011-972um-1224um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161014-1296um-1476um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161015-1368um-1476um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161016-1368um-1476um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161018-1552um-2304um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161019-1624um-2340um-161019-095311_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161020-1768um-2376um-161020-102642_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161021-1840um-2412um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161023-1872um-2448um-161023-123504_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161025-1944um-2232um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161026-1944um-2232um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161027-1944um-2232um-161027-194125_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161030-2052um-2232um-161030-111142_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161031-2088um-2232um-161031-104001_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161101-2088um-2232um-merge_behavior+ecephys.nwb",
 "sub-DT5/sub-DT5_ses-20161104-2124um-2232um-merge_behavior+ecephys.nwb",
],
"DT7": [
 "sub-DT7/sub-DT7_ses-20170315-288um-0um-170315-191240_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170317-288um-0um-merge_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170318-288um-0um-170318-115923_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170320-432um-0um-merge_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170321-216um-0um-170321-120259_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170324-576um-144um-170324-123932_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170325-576um-360um-170325-125511_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170401-864um-936um-merge_behavior+ecephys.nwb",
 "sub-DT7/sub-DT7_ses-20170412-1440um-1170um-170412-103153_behavior+ecephys.nwb",
],
}
