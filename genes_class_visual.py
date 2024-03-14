#!/usr/bin/env python
# coding: utf-8


import os
from collections import Counter
import pandas as pd
import numpy as np
import chart_studio.plotly as py
import cufflinks as cf
import seaborn as sns
import plotly.express as px
import plotly.graph_objs as go
import regex as re
from collections import defaultdict
get_ipython().run_line_magic('matplotlib', 'inline')

# Make Plotly work in your Jupyter Notebook
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode(connected=True)
# Use Plotly locally
cf.go_online()

# Путь к директории с файлами кластеров
cluster_dir = "../../common/clastering_hole_exempls/clusters/"
# Путь к файлу с координатами генов
genes_path = "../../../avel/common/clastering_hole_exempls/genes_in_all_exempl.csv"

unique_samples = ['drad_1778', 'drad_200', 'drad_215', 'duni_224', 'duni_242', 'duni_243', 'dval_218', 'dval_245']

# Функция для подсчета количества вхождений образца в кластеры
def categorize_sample(sample, cluster_df):
    # Подсчёт вхождений образца в каждый кластер
    sample_counts = Counter()
    for _, row in cluster_df.iterrows():
        sample_counts[row['cluster_name']] += row['members'].count(sample)

    # Категоризация образца по вхождениям в кластеры
    if sample_counts['rudis'] == 1 and sample_counts['caucasica'] == 0:
        return 'R'
    elif sample_counts['rudis'] == 2 and sample_counts['caucasica'] == 0:
        return 'RR'
    elif sample_counts['caucasica'] == 1 and sample_counts['rudis'] == 0:
        return 'C'
    elif sample_counts['caucasica'] == 2 and sample_counts['rudis'] == 0:
        return 'CC'
    elif sample_counts['caucasica'] == 1 and sample_counts['rudis'] == 1:
        return 'RC'
    elif sample_counts['undefind'] == 1 and sample_counts['rudis'] == 1:
        return 'RN'
    elif sample_counts['caucasica'] == 1 and sample_counts['undefind'] == 1:
        return 'CN'
    elif sample_counts['undefind'] == 1:
        return 'N'
    elif sample_counts['undefind'] >= 2:
        return 'NN'
    elif sample_counts['caucasica'] == 0 and sample_counts['undefind'] == 0 and sample_counts['rudis'] == 0:
        return 'O'
    else:
        return 'X'


# Считывание датафрейма с координатами генов
genes_df = pd.read_csv(genes_path)

# Проход по всем файлам кластеров и категоризация образцов
for file_name in os.listdir(cluster_dir):
    if file_name.endswith("_clusters.csv"):
        gene_id = file_name.split('_')[0]
        cluster_df = pd.read_csv(os.path.join(cluster_dir, file_name))

        # Категоризация каждого уникального образца и обновление датафрейма координат
        for sample in unique_samples:
            category = categorize_sample(sample, cluster_df)
            # Предполагаем, что в genes_df есть столбец с именем образца для обновления категории
            genes_df.loc[genes_df['busco_id'] == gene_id, sample] = category

# Сохранение обновленного датафрейма в новый CSV файл
# genes_df.to_csv('categorized_genes.csv', index=False)

# Define the data as a list of lists
data = [
    ["NC_046312.1", 1, 133750839, "1"],
    ["NC_046313.1", 1, 120754861, "2"],
    ["NC_046314.1", 1, 118060892, "3"],
    ["NC_046315.1", 1, 100704004, "4"],
    ["NC_046316.1", 1, 95883647, "5"],
    ["NC_046317.1", 1, 95499011, "6"],
    ["NC_046318.1", 1, 86565987, "7"],
    ["NC_046319.1", 1, 84102018, "8"],
    ["NC_046320.1", 1, 74758939, "9"],
    ["NC_046321.1", 1, 70688303, "10"],
    ["NC_046322.1", 1, 61599656, "11"],
    ["NC_046323.1", 1, 56660093, "12"],
    ["NC_046324.1", 1, 52485139, "13"],
    ["NC_046325.1", 1, 51720411, "14"],
    ["NC_046326.1", 1, 43910048, "15"],
    ["NC_046327.1", 1, 40449263, "16"],
    ["NC_046328.1", 1, 40056199, "17"],
    ["NC_046329.1", 1, 12095641, "18"],
    ["NC_046331.1", 1, 47440541, "Z"],
    ["NC_046330.1", 1, 3982059, "W"]
]

# Create a DataFrame
scaffold_df = pd.DataFrame(data, columns=["scaffold", "start", "end", "chr"])

fig = go.Figure()
chr_width = 0.1
mark_width = 0.2
sample = 'duni_224'

# Функция для преобразования значений в 'chr' из строк в числа с учетом 'Z' и 'W'
def convert_chr_to_numeric(chr_value):
    if chr_value == 'Z':
        return 19  # Например, присвоить Z значение 19
    elif chr_value == 'W':
        return 20  # Например, присвоить W значение 20
    else:
        return int(chr_value)  # Остальные значения преобразуются в числа

# Применяем функцию к столбцу и добавляем смещение
genes_df['chr_numeric'] = genes_df['chr'].apply(convert_chr_to_numeric)
scaffold_df['chr_numeric'] = scaffold_df['chr'].apply(convert_chr_to_numeric)

# Create scatter trace of Chromosomes
fig.add_trace(go.Bar(
    x=scaffold_df['end'],
    y=scaffold_df['chr_numeric'],
    orientation='h',
    name = 'Chromosome',
    width = chr_width,
))

fig.update_layout(barmode='overlay')

# Title of figure and axises
fig.update_layout(title={
        'text': "Gene types",
        'y':0.99,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'}, 
        xaxis_title='Bp', yaxis_title='Chromosome name')


fig.update_layout(
    # Shows gray line without grid, styling fonts, linewidths and more
    xaxis=dict(
        showline=True,
        showgrid=False,
        showticklabels=True,
        linecolor='rgb(204, 204, 204)',
        linewidth=1,
        ticks='outside',
        rangemode="nonnegative",
        tickfont=dict(
            family='Arial',
            size=15,
            color='rgb(82, 82, 82)',
        ),
    ),
    # Turn off everything on y axis
    yaxis=dict(
        showgrid=False,
        zeroline=False,
        showline=False,
        showticklabels=True,
        ticklabelstep=1,
        tickwidth=15,
        tickfont=dict(
            family='Arial',
            size=15,
            color='rgb(82, 82, 82)',
        ),
    ),
    width=1400,
    height=900,
    margin=dict(
        autoexpand=True,
        #l=150,
        #r=20,
        #t=110,
    ),
    showlegend=True,
    plot_bgcolor='white'
)

fig.update_layout(legend = dict(font = dict(family = "Arial", size = 15, color = "black")))


fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='RC']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='RC']['length']*3,
    y=genes_df.loc[genes_df[sample]=='RC']['chr_numeric'],  
    marker_color="Black",
    orientation="h",
    width=mark_width,
    name="RC - Rudis and Caucasica"
))

fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='RR']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='RR']['length']*3,
    y=genes_df.loc[genes_df[sample]=='RR']['chr_numeric'] + 0.2,  
    marker_color="Red",
    orientation="h",
    width=mark_width,
    name="RR - Double Rudis"
))

fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='R']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='R']['length']*3,
    y=genes_df.loc[genes_df[sample]=='R']['chr_numeric'] + 0.4, 
    marker_color="Maroon",
    orientation="h",
    width=mark_width,
    name="R - Only one Rudis"
))

fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='CC']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='CC']['length']*3,
    y=genes_df.loc[genes_df[sample]=='CC']['chr_numeric'] - 0.2,  
    marker_color="Lime",
    orientation="h",
    width=mark_width,
    name="CC - Double Caucasica"
))

fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='C']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='C']['length']*3,
    y=genes_df.loc[genes_df[sample]=='C']['chr_numeric'] - 0.4,  
    marker_color="Green",
    orientation="h",
    width=mark_width,
    name="C - Only one Caucasica"
))

fig.add_trace(go.Bar(
    base=genes_df.loc[genes_df[sample]=='O']['gene_start'],
    x=genes_df.loc[genes_df[sample]=='O']['length']*3,
    y=genes_df.loc[genes_df[sample]=='O']['chr_numeric'],  
    marker_color="White",
    orientation="h",
    width=0.1,
    name="O - Losed"
))

fig.show()

