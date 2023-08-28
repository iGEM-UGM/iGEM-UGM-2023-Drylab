import RNA
import pandas as pd

def mfe_segment(data):
  list_domain = ['ABS',
                  'ABS_',
                  'ATS',
                  'Complementstart',
                  'DBSCom',
                  'Loop',
                  'RBS',
                  'StartCodon',
                  '_DBS',
                  'loopCom',
                  'loopRBS',
                  'miR21bottomseal1',
                  'miR21bottomseal2',
                  'miR21topseal',
                  'miR92bottomseal1',
                  'miR92bottomseal2',
                  'miR92topseal',
                  'toeholdmiR92',
                  'toeholdmir21']
  list_segment = [['ABS', 'Complementstart', 'ABS_', 'Loop', 'ATS', 'loopRBS', 'miR21bottomseal1'],
                  ['toeholdmir21', 'miR21topseal', 'miR21bottomseal2', 'miR92bottomseal1', 'toeholdmiR92', 'miR92topseal'],
                  ['miR92bottomseal2', 'RBS', 'loopCom', '_DBS', 'StartCodon', 'DBSCom']]
  res = dict()
  list_name = ['ABS-miR21bottomseal1', 'toeholdmir21-miR92topseal', 'miR92bottomseal2-end']
  for j, name in zip(list_segment, list_name):
    # get The RNA sequence segment
    seq = ''
    for i in j:
      # print(i)
      dom = data.loc[data['Domain_Strand'] == str(i), 'Sequence'].values[0]
      seq += dom
    seq = seq.lower()
    print(seq)
    # compute minimum free energy (MFE) and corresponding structure
    (ss, mfe) = RNA.fold(seq)
    res['mfe_' + name] = mfe
  return(res)

def main():
  table_path = input("Input your table get_sequence output: " )
  df = pd.read_csv('./'+table_path)
  datas = []
  for i in df['Trial'].unique().tolist():
      datas.append(mfe_segment(df[df['Trial']==i]))

if __name__ == "__main__":
  main()