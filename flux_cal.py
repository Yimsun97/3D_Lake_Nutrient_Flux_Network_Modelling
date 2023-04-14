import re
import numpy as np
import pandas as pd

from meta import *

scenarios = range(11)  # ['S0','S1','S2','S3','S3a']

atm_deso = pd.read_excel(f'{input_data_path}/ATM_deso.xlsx')
with open(f'{input_data_path}/wq3dwc.inp', 'r', encoding='utf-8') as r:
    wq3dwc_lines = r.readlines()

all_fluxes = ['JDAY', 'ALGAEN', 'ALGAEN_fix', 'ALGAEN_SET1', 'ALGAEN_PS', 'ALGAEN_EXIT',
              'ALGAEP', 'ALGAEP_SET1', 'ALGAEP_PS', 'ALGAEP_EXIT',
              'PO4', 'PO4_BEN1', 'PO4_SET1', 'PO4_PS', 'PO4_EXIT', 'PO4_ALGAE_UPTK', 'PO4_FROM_DOP',
              'PO4_FROM_ALGAE', 'PO4_RESUS',
              'DOP', 'DOP_PS', 'DOP_EXIT', 'DOP_FROM_ALGAE', 'DOP_FROM_POP',
              'POP', 'POP_SET1', 'POP_PS', 'POP_EXIT', 'POP_FROM_ALGAE', 'POP_RESUS',
              'NH4', 'RNH4_BEN1', 'RNH4_PS', 'RNH4_EXIT', 'RNH4_ALGAE_UPTK', 'NH4_FROM_DON',
              'NH4_FROM_ALGAE', 'NH4_RESUS',
              'NO3', 'NO3_denit', 'RNO3_BEN1', 'RNO3_PS', 'RNO3_EXIT', 'RNO3_ALGAE_UPTK',
              'NO3_FROM_NH4', 'NO3_RESUS',
              'DON', 'DON_PS', 'DON_EXIT', 'DON_FROM_ALGAE', 'DON_FROM_PON',
              'PON', 'PON_SET1', 'PON_PS', 'PON_EXIT', 'PON_FROM_ALGAE', 'PON_RESUS',
              'N_atm', 'N_resus', 'P_atm', 'P_resus', 'SED_P', 'SED_N']
export_fluxes = ['JDAY', 'ALGAEN', 'ALGAEP', 'PO4', 'DOP', 'POP', 'NH4', 'NO3', 'DON', 'PON',
                 'SED_P', 'SED_N']

# get parameters of algal dynamics
for i, line in enumerate(wq3dwc_lines):
    line = line.strip()
    if line[0:3] == 'C08':
        lineN0 = i
    if line[0:3] == 'C34':
        lineN1 = i
        break

wq3dwc = wq3dwc_lines[lineN0: lineN1]
para = {}
for i, line in enumerate(wq3dwc):
    line = line.strip()
    # parameter lines if starting with digits
    if re.search('^[0-9\.]', line):
        line_name = re.split('[ \s\t]{1,}', wq3dwc[i - 1].strip())[1:]
        line_value = re.split('[ \s\t]{1,}', line)
        for k, val in enumerate(line_value):
            para[line_name[k].upper()] = float(val)

for scn in scenarios:
    N = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQ_NBALANCE.OUT', sep='\s+')
    P = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQ_PBALANCE1.OUT', sep='\s+')
    KINE = pd.read_table(f'{scn_path}\\S{scn:02d}\\WQ_KINE_FLUX.OUT', sep='\s+')
    N = N.rename(columns={'NFIXING': 'ALGAEN_fix', 'DENITRIFICATION': 'NO3_denit', 'ATM_DEPO': 'N_atm'})
    P = P.rename(columns={'TPOP_MASS': 'POP', 'TDOP_MASS': 'DOP', 'ATM_DEPO': 'P_atm'})
    KINE = KINE.rename(columns={'LDP_FROM_ALGAE': 'DOP_FROM_ALGAE', 'P_ALGAE_UPTK': 'PO4_ALGAE_UPTK',
                                'DOP_FROM_LOP': 'DOP_FROM_POP'})

    flux_inp = pd.concat([N, P.iloc[:, 1:], KINE.iloc[:, 1:]], axis=1)
    names = [i.replace('_MASS', '') for i in flux_inp.columns]
    new_name = dict(zip(flux_inp.columns, names))
    flux_inp = flux_inp.rename(columns=new_name)
    flux_inp['PON_RESUS'] = flux_inp['LON_RESUS'] + flux_inp['RON_RESUS']
    flux_inp['POP_RESUS'] = flux_inp['LOP_RESUS'] + flux_inp['ROP_RESUS']

    if scn == 10:
        sed_n_init = 70124.33024
        sed_p_init = 28372.34006

    # add stock 'SED_N' and 'SED_P'
    flux_inp.loc[:, 'SED_N'] = sed_n_init - flux_inp['ALGAEN_SET1'] - flux_inp['PON_SET1'] - \
                               flux_inp['RNH4_BEN1'] - flux_inp['RNO3_BEN1'] - \
                               flux_inp['NH4_RESUS'] - flux_inp['PON_RESUS'] - flux_inp['NO3_RESUS']

    flux_inp.loc[:, 'SED_P'] = sed_p_init - flux_inp['ALGAEP_SET1'] - flux_inp['POP_SET1'] - \
                               flux_inp['PO4_BEN1'] - flux_inp['PO4_RESUS'] - flux_inp['POP_RESUS']

    flux_inp = flux_inp.loc[:, all_fluxes]

    flux_inp = flux_inp.loc[flux_inp['JDAY'] % 1 < 0.002].copy()
    flux_inp.loc[:, 'JDAY'] = flux_inp.JDAY.astype(int)
    flux_inp = flux_inp.set_index(flux_inp.loc[:, 'JDAY']).copy()

    # update NH4_FROM_ALGAE based on PO4 from algae
    ALGAE_BM = (flux_inp['PO4_FROM_ALGAE'] * para['CPPRM1'] /
                (para['BMRC'] * para['FPIC'] + para['PRRC'] * para['FPIP']) * para['BMRC'])
    ALGAE_PR = (flux_inp['PO4_FROM_ALGAE'] * para['CPPRM1'] /
                (para['BMRC'] * para['FPIC'] + para['PRRC'] * para['FPIP']) * para['PRRC'])
    NH4_ALGAE_BM = ALGAE_BM * para['ANCC'] * para['FNIC']
    NH4_ALGAE_PR = ALGAE_PR * para['ANCC'] * para['FNIP']
    flux_inp.loc[:, 'NH4_FROM_ALGAE'] = NH4_ALGAE_BM + NH4_ALGAE_PR

    # calculation daily difference
    # initial values for the 1st day (JDAY=2193)
    # difference values for the days after the 1st day
    JDAY = flux_inp['JDAY']
    mass = flux_inp.loc[:, export_fluxes]
    data1 = flux_inp.copy()
    data1.index = range(flux_inp.index[0] + 1, flux_inp.index[-1] + 2, 1)
    data2 = flux_inp - data1
    flux_inp.update(data2)
    flux_inp.loc[:, 'JDAY'] = JDAY

    # add atmospheric deposition
    atm_deso.index = atm_deso['TIME']
    flux_inp['PON_ATM'] = atm_deso['RPON_ATM'] + atm_deso['LPON_ATM']
    flux_inp['DON_ATM'] = atm_deso['RDON_ATM'] + atm_deso['LDON_ATM']
    flux_inp['NH4_ATM'] = atm_deso['NH4_ATM']
    flux_inp['NO3_ATM'] = atm_deso['NO3_ATM']
    flux_inp['POP_ATM'] = atm_deso['RPOP_ATM'] + atm_deso['LPOP_ATM']
    flux_inp['DOP_ATM'] = atm_deso['RDOP_ATM'] + atm_deso['LDOP_ATM']
    flux_inp['PO4_ATM'] = atm_deso['PO4_ATM']

    # aggregation statistics
    flux_day = flux_inp.copy()
    flux_day['TIME'] = pd.to_datetime(jday_start) + flux_inp['JDAY'].apply(lambda x: pd.Timedelta(days=x-1))
    flux_day['year'] = flux_day['TIME'].dt.year
    flux_day['month'] = flux_day['TIME'].dt.month
    flux_day['week'] = flux_day['TIME'].dt.isocalendar().week
    flux_day['day'] = flux_day['TIME'].dt.day
    flux_day['yr_mon'] = flux_day['year'].map(str) + '_' + flux_day['month'].apply(lambda x: '%02d' % x)
    flux_day['yr_week'] = flux_day['year'].map(str) + '_' + flux_day['week'].apply(lambda x: '%02d' % x)

    columns = ['TIME', 'year', 'month', 'week', 'day', 'yr_mon', 'yr_week']
    mass.loc[:, columns] = flux_day.loc[:, columns].values
    mass.to_excel(f'{flux_output_path}\\mass_NP_S{scn:02d}.xlsx', index=False)
    flux_day.drop(['ALGAEN', 'ALGAEP', 'PO4', 'DOP', 'POP', 'NH4', 'NO3', 'DON', 'PON'], axis=1).to_excel(
        f'{flux_output_path}\\flux_NP_S{scn:02d}.xlsx', index=False)
    print(scn)
