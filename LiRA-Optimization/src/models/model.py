def get_df(name):
  url = f'https://raw.githubusercontent.com/iGEM-UGM/iGEM-UGM-drylab/main/{name}'
  response = requests.get(url)
  with open(name, 'wb') as f:
      f.write(response.content)

  return(pd.read_csv(name))

import pandas as pd
import numpy as np

df1 = get_df('nupack_10_weighted.csv')
df2 = get_df('mfe_data_10trial_revisi.csv')
df = pd.concat([df1, df2], axis=1)

# get columns that has variation of complexes
cols_complex = []
for i in df.columns[df.columns.str.startswith('on_level_')].tolist()[1:-1]:
  cols_complex.append(i.split('on_level_')[1])

# get the list of columns to get into regression
  _param = ['_mfe',
          '_pfunc',
          '_free_energy',
          '_ensemble_size',
          '_mfe_energy',
          '_mfe_stack_energy',
          '_t1',
          '_t1_rank']
  param_  = ['on_level_',
            'off_level_',
            'on_mean_',
            'on_sum_',
            'on_stdev_',
            'off_mean_',
            'off_sum_',
            'off_stdev_']
  param = ['mfe_ABS-miR21bottomseal1', 'mfe_toeholdmir21-miR92topseal',
          'mfe_miR92bottomseal2-end',
          'gc_content', 'structure_prob', 'defect']

def prepro(df, mode = 'infinity'):
  # change datatype to float
  df['(a)_pfunc'] = df['(a)_pfunc'].astype(float)
  df['(a)_ensemble_size'] = df['(a)_ensemble_size'].astype(float)
  df['(a+c+b)_ensemble_size'] = df['(a+c+b)_ensemble_size'].astype(float)
  df['(a+b+c)_ensemble_size'] = df['(a+b+c)_ensemble_size'].astype(float)
  df['(a+b)_ensemble_size'] = df['(a+b)_ensemble_size'].astype(float)
  df['(a+c)_ensemble_size'] = df['(a+c)_ensemble_size'].astype(float)
  df['(a+a)_ensemble_size'] = df['(a+a)_ensemble_size'].astype(float)
  df_non_inf = pd.DataFrame()
  if mode == 'infinity':
    # Check which rows contain infinity values
    rows_with_inf = df.isin([np.inf]).any(axis=1)

    # Drop rows with infinity values
    df_non_inf = df.drop(index = pd.Series(rows_with_inf)[pd.Series(rows_with_inf)].index[0])
  return(df, df_non_inf)

def get_param(df, _param, cols_complex):

  params_all = []
  for i in cols_complex:
    for j in _param:
      params_all.append(str(i+j))
    # for k in param_:
    #   params_all.append(str(k+i))
  params_all.extend(param)
  return(params_all)
print(len(cols_complex)*(len(_param))+len(param))
len(params_all)

# modelling
def get_mse_rmse(y_test, y_pred):
  mse = mean_squared_error(y_test, y_pred)
  rmse = np.sqrt(mse)
  return mse, rmse

def linreg(X, y):
  sc_X = StandardScaler()
  X = sc_X.fit_transform(X)
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

  # Build linear regression model
  model = LinearRegression()
  model.fit(X_train, y_train)
  # Make predictions on testing set
  y_pred = model.predict(X_test)

  mse, rmse = get_mse_rmse(y_test, y_pred)

  print('MSE:', mse)
  print('RMSE:', rmse)
  return(model)

# Import necessary libraries
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

# Split data into training and testing sets
X = df_non_inf[params_all]
y = df_non_inf['on_off_ratio_avg']
model = linreg(X, y)

def display_metrics(sf, sc):
  display(pd.DataFrame({
      'Features': sf,
      'Coefficients': sc,
      'Absolute Coefficients': np.abs(sc)
  }))

def metrics(model, params_all):
  # Get the coefficients of the linear regression model
  coefficients = model.coef_

  # Sort the coefficients in descending order
  sorted_indices = np.argsort(coefficients)[::-1].tolist()
  print(sorted_indices)
  sorted_coefficients = coefficients[sorted_indices]
  sorted_features = sorted(params_all, key=lambda x: sorted_indices.index(params_all.index(x)))

  # Print or visualize the feature importances
  for feature, importance in zip(sorted_features, sorted_coefficients):
      print(f"{feature}: {importance}")

  # Create a bar plot of feature importances
  plt.bar(range(len(sorted_features)), sorted_coefficients)
  plt.xticks(range(len(sorted_features)), sorted_features, rotation='vertical')
  plt.xlabel('Features')
  plt.ylabel('Coefficient')
  plt.title('Feature Coefficients')
  plt.show()
  return(zip(sorted_features, sorted_coefficients))
