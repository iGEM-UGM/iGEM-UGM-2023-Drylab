import nupack
from nupack import *
def get_rbs_index(c):
    a = ''
    for i in c.strands:
        a += str(i)
    rbs = "UAGAGGAGAUG"
    index = a.find(rbs)
    return(int(index))

def calculate_gc_percentage(string):
    gc_count = string.count("G") + string.count("C")
    total_count = len(string)
    gc_percentage = (gc_count / total_count) * 100

    return(gc_percentage)

import re
import numpy as np
def check_illegal_sites(sequence, enzymes):
    illegal_sites = 0
    for enzyme, site in enzymes.items():
        matches = re.finditer(site, sequence)
        for match in matches:
            illegal_sites += 1

    return(illegal_sites)

def get_index(text):
    periods = []
    parentheses = []

    for match in re.finditer(r'\.|[()]', text):
        if match.group() == '.':
            periods.append(match.start())
        elif match.group() == '(' or match.group() == ')':
            parentheses.append(match.start())
    return periods, parentheses

def get_prob(target, matrix):
    # target = str(tube_results.complexes[Complex([b, a, c], name='(b+a+c)')].mfe[0][0]).replace("+", "")
    target = target.replace('+', '')
    # matrix = pairs(Complex([b, a, c], name='(b+a+c)'), model=model1).to_array()
    diagonal_elements = np.diag(matrix)
    non_diagonal_array = matrix - diagonal_elements
    prob = []
    for i in range(len(target)):
        if target[i] == ')' or target[i] == '(':
            prob.append(np.max(non_diagonal_array[i]))
        elif target[i] == '.':
            prob.append(diagonal_elements[i])
        else:
            continue
    return prob

def get_probability_onoff(index, data):
    prob = []
    for i in index:
        prob.append(data[i])
    return prob

def on_off(text, data):
    periods, parentheses = get_index(text)
    prob_on  = []
    prob_off = []
    prob_on = get_probability_onoff(periods, data)
    prob_off = get_probability_onoff(parentheses, data)
    return(prob_on, prob_off)

def on_off_level(text, data):
    prob_on, prob_off = on_off(text, data)
    on_level = sum(prob_on)/len(text)
    off_level = sum(prob_off)/len(text)
    return(on_level, off_level)

def stats(lists):
    try:
        mean = sum(lists) / len(lists)
        squared_differences = [(x - mean)**2 for x in lists]
        variance = sum(squared_differences) / len(lists)
        standard_deviation = variance**0.5
        sum_res = sum(lists)

        return(mean, sum_res, standard_deviation)
    except ZeroDivisionError:
        print("Error: Division by zero.")
        return 0, 0, 0

def get_analyze(t, df):
    seq = str.upper((df.loc[(df['Trial'] == t) & (df['Domain_Strand'] == 'LiRA'), 'Sequence'].values)[0])
    mir21 = 'UAGCUUAUCAGACUGAUGUUGA'
    mir92a = 'UAUUGCACUUGUCCCGGCCUGU'


    # define illegal sites enzymes
    '''
    Ecori: GAATTC GAAUUC
    Xbai: TCTAGA UCUAGA
    Spei: ACTAGT ACUAGU
    Psti: CTGCAG CUGCAG
    noti:GCGGCCGC
    '''
    enzymes = {
      "EcoRI": "GAATTC",
      "EcoRI2": "GAAUUC",
      "XbaI": "TCTAGA",
      "XbaI2": "UCUAGA",
      "SpeI": "ACTAGT",
      "SpeI2": "ACUAGU",
      "PstI": "CTGCAG",
      "PstI2": "CUGCAG",
      "NotI": "GCGGCCGC"
    }
    # analysis job

    # specify strands
    a = Strand(seq, name='a')
    b = Strand(mir21, name='b')
    c = Strand(mir92a, name='c')

    # specify tubes
    t1 = Tube(strands={a: 1e-6, b: 1e-6, c: 1e-6}, complexes=SetSpec(max_size=3, include=[[a,b], [a,c], [a,b,c]]), name='t1')
    # t2 = Tube(strands={a: 1e-10, b: 1e-9}, complexes=SetSpec(max_size=2), name='t2')

    # analyze tubes
    model1 = Model()
    tube_results = tube_analysis(tubes=[t1], model=model1, compute=['mfe','ensemble_size'])

    # Split the ASCII result into lines
    lines = str(tube_results).split('\n')

    # Extract the column names
    columns = lines[0].split()

    # Create an empty list to store the data rows
    concen_data = []

    # Iterate over the lines starting from the second line
    for line in lines[1:]:
        # Split the line into individual values
        values = line.split()
        # Append the values as a row to the data list
        concen_data.append(values)

    # Iterate through the list and find the index containing 'Concentration'
    index = None
    for i, sublist in enumerate(concen_data):
        if 'Concentration' in sublist:
            index = i
            break

    new_data = concen_data[index+1:]
    # Create the dataframe
    concen = pd.DataFrame(new_data[1:], columns=new_data[0][:-1], )

    #  data get from each sequence
    c0 = Complex([a], name='(a)')
    c1 = Complex([a, c, b], name='(a+c+b)')
    c2 = Complex([a, b, c], name='(a+b+c)')
    c3 = Complex([a, b], name='(a+b)')
    c4 = Complex([a, c], name='(a+c)')
    c5 = Complex([a, a], name='(a+a)')

    data = dict ()
    data['sequence'] = seq
    for c in [c0, c1, c2, c3, c4, c5]:
        data[str(c.name) + '_mfe'] =  tube_results.complexes[c].mfe_stack
        data[str(c.name) + '_target']  = str(tube_results.complexes[c].mfe[0][0])
        data[str(c.name) + '_pfunc'] =  tube_results.complexes[c].pfunc
        data[str(c.name) + '_free_energy'] =  tube_results.complexes[c].free_energy
        data[str(c.name) + '_ensemble_size'] =  tube_results.complexes[c].ensemble_size
        data[str(c.name) + '_mfe_energy'] =  tube_results.complexes[c].mfe_stack
        data[str(c.name) + '_mfe_stack_energy'] =  tube_results.complexes[c].mfe[0][2]
        data[str(c.name) + '_t1'] = concen.loc[concen['Complex'] == str(c.name), 't1'].values[0]
        data[str(c.name) + '_t1_rank'] = concen.index[concen['Complex'] == str(c.name)][0]
        data['on_' + str(c.name)] = len(re.findall(r'\.', str(tube_results.complexes[c].mfe[0][0]))[get_rbs_index(c):get_rbs_index(c)+18])
        print(get_rbs_index(c))
        data['off_'+ str(c.name)] = len(re.findall(r'[()]', str(tube_results.complexes[c].mfe[0][0]))[get_rbs_index(c):get_rbs_index(c)+18])
        data['prob_matrix_'+ str(c.name)] = pairs(c, model=model1).to_array()

        data['prob_' + str(c.name)] = get_prob(str(tube_results.complexes[c].mfe[0][0]).replace("+", ""), pairs(c, model=model1).to_array())
        data['prob_on_' + str(c.name)], data['prob_off_' + str(c.name)] = on_off(str(tube_results.complexes[c].mfe[0][0])[get_rbs_index(c):get_rbs_index(c)+18], data['prob_' + str(c.name)])
        data['on_level_' + str(c.name)], data['off_level_' + str(c.name)] = on_off_level(str(tube_results.complexes[c].mfe[0][0])[get_rbs_index(c):get_rbs_index(c)+18], data['prob_' + str(c.name)])
        data['on_mean_' + str(c.name)], data['on_sum_' + str(c.name)], data['on_stdev_' + str(c.name)] = stats(data['prob_on_' + str(c.name)])
        data['off_mean_' + str(c.name)], data['off_sum_' + str(c.name)], data['off_stdev_' + str(c.name)] = stats(data['prob_off_' + str(c.name)])

    data['defect'] = defect(strands=[seq], structure=data['(a)_target'], model=model1)
    data['structure_prob'] = structure_probability(strands=[seq], structure=data['(a)_target'], model=model1)
    # data['prob_matrix'] = pairs(strands=[seq], model=model1).to_array()
    # data['prob_diagonal'] = np.diagonal(pairs(strands=[seq], model=model1).to_array())
    data['gc_content'] = calculate_gc_percentage(seq)
    data['illegal_count'] = check_illegal_sites(seq, enzymes)
    return(data)

def get_on_off_each_seq(cols_complex, result)
  for i in cols_complex:
      print(i)
      on_off_ratio = []
      on_off_minus = []
      for j in range(len(result)):
          # print(result['on_level_' + str(i)][j])
          # print(result['off_level_' + str(i)][j])
          ratio = result['on_level_' + str(i)][j] / result['off_level_' + str(i)][j]
          minus = result['on_level_' + str(i)][j] - result['off_level_' + str(i)][j]
          # print(ratio)
          on_off_ratio.append(ratio)
          on_off_minus.append(minus)
  return(on_off_ratio, on_off_minus)

def get_metrics(cols_complex, result):
  ratio_avg_list = []
  minus_avg_list = []
  on_level_avg_list = []
  off_level_avg_list = []
  w_ratio_avg_list = []
  w_minus_avg_list = []
  w_on_level_avg_list = []
  w_off_level_avg_list = []
  for j in range(len(result)):
      ratio_avg = 0
      minus_avg = 0
      on_level_avg = 0
      off_level_avg = 0
      weighted_on_level_avg = 0
      weighted_off_level_avg = 0
      weighted_minus_avg = 0
      weighted_ratio_avg = 0
      for i in cols_complex:
          # complex regression: ['(a+c+b)', '(a+b+c)', '(a+b)', '(a+c)']
          ratio_avg += result['on_off_ratio_' + str(i)][j]
          minus_avg += result['on_off_minus_' + str(i)][j]
          on_level_avg += result['on_level_' + str(i)][j]
          off_level_avg += result['off_level_' + str(i)][j]
          weighted_on_level_avg += (1/(result[str(i) + '_t1_rank'][j]+1))*result['on_level_' + str(i)][j]
          weighted_off_level_avg += (1/(result[str(i) + '_t1_rank'][j]+1))*result['off_level_' + str(i)][j]
          weighted_ratio_avg += (1/(result[str(i) + '_t1_rank'][j]+1))*result['on_off_ratio_' + str(i)][j]
          weighted_minus_avg += (1/(result[str(i) + '_t1_rank'][j]+1))*result['on_off_minus_' + str(i)][j]
          # print(result['on_off_ratio_' + str(i)][j])
      ratio_avg_list.append(ratio_avg/len(cols_complex))
      minus_avg_list.append(ratio_avg/len(cols_complex))
      on_level_avg_list.append(on_level_avg/len(cols_complex))
      off_level_avg_list.append(off_level_avg/len(cols_complex))
      #weighted
      w_minus_avg_list.append(weighted_ratio_avg/len(cols_complex))
      w_ratio_avg_list.append(weighted_ratio_avg/len(cols_complex))
      w_on_level_avg_list.append(weighted_on_level_avg/len(cols_complex))
      w_off_level_avg_list.append(weighted_off_level_avg/len(cols_complex))
      return(ratio_avg_list, minus_avg_list, on_level_avg_list, off_level_avg_list,
             w_ratio_avg_list, w_minus_avg_list, w_on_level_avg_list, w_off_level_avg_list)
  
def get_complex(df):
  cols_complex = []
  for i in df.columns[df.columns.str.startswith('on_level_')].tolist()[1:-1]:
    cols_complex.append(i.split('on_level_')[1])
  print(cols_complex)
  return(cols_complex)

import pandas as pd
def main():
    table_path = input("Input your table get_sequence output: " )
    df = pd.read_csv('./'+table_path)
    datas = []
    for i in df['Trial'].unique().tolist():
        datas.append(get_analyze(i, df))
    result = pd.DataFrame(datas)

    cols_complex = get_complex(result)
    result['on_off_ratio_' + str(i)], result['on_off_minus_' + str(i)] = get_on_off_each_seq(cols_complex, result)
    result['on_off_ratio_avg'], result['on_off_minus_avg'], result['on_level_avg'], result['off_level_avg'], result['weighted_on_off_ratio_avg'], result['weighted_on_off_minus_avg'], result['weighted_on_level_avg'], result['weighted_off_level_avg'] = get_metrics(cols_complex, result)
    result.to_csv('nupack_params.csv')

if __name__ == "__main__":
  main()