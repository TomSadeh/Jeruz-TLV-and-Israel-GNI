# -*- coding: utf-8 -*-

# Importing the required libraries.

import pandas as pd
import pyreadstat as prs
import numpy as np
import matplotlib.pyplot as plt

# Defining a small and simple function to invert strings for the plots and better readablity.

def invert(string):
    """
    A function which invert a string.
    Parameters
    ----------
    string : string
        the string to invert.

    Returns
    -------
    string
        An inverted string.

    Required libraries
    ------------------
    None.
    """
    return string[::-1]

# Importing a file with the names of the household files, and a file with CPI data. Enter the files location here.

file_names = pd.read_csv('new_file_names.csv', index_col = 'year')  
cpi = pd.read_csv('cpi_for_git.csv', index_col = 'year')

# Creating a small dictionary with the period of the analysis.

period = {'start' : 1997,
          'end' : 2019}

#Creating a variable to contain the base address of the files. Enter yours here.

base_add = r''

# Creating an empty DataFrame to contain the results of the analysis.

results = pd.DataFrame()

# The main loop of the analysis. In the loop we load each year files and make a calculation each year and saving it to the results DataFrame.

for year, folder, file_mbs, file_prod in zip(np.arange(period['start'], period['end'] + 1),
                                                        file_names.loc[period['start']:period['end'], 'folder_address'],
                                                        file_names.loc[period['start']:period['end'], 'mb'],
                                                        file_names.loc[period['start']:period['end'], 'prod']):
    
    # 1997-2003 are years with differnet file type, so we load them accordingly with pyreadstat.
    
    if year < 2004:
        mbs = prs.read_dta(base_add + '\\' + folder + '\\' + file_mbs + '.dta', encoding = 'windows-1255')
        prod = prs.read_dta(base_add + '\\' + folder + '\\' + file_prod + '.dta', encoding = 'windows-1255')
        
        # Taking the column headers from the metadata, as the defult ones aren't readable.
        
        mbs[0].columns = mbs[1].column_labels
        prod[0].columns = prod[1].column_labels
        
        # Renaming some of the hebrew headed columns to match the newer files, and setting the household number as the index for the households DataFrame for merging with the product DataFrame later.
        
        mbs[0].rename(columns = {'מספר משפח' : 'misparmb', 'משקל משפח' : 'weight', 'גודל משפח' : 'nefashot', 'צורת ישו' : 'zurat_y'}, inplace = True)
        
        # 1997-2002 contains a pivoted prod file. 2003 does not, so we pivot the DataFrame in this year. 
        
        if year < 2003:
            prod[0].set_index('misparmb', inplace = True)
            prod_pivoted = prod[0]
                
        else:
            prod_pivoted = prod[0].pivot(index = 'misparmb', columns = 'prodcode', values = 'schum')

        # After we finished with the metadata (it's stored in mbs[1]), we discard it, for easeier accessing the DataFrame later.
    
        mbs = mbs[0]
        
        # 2004 is in a differnet format that the others, so we need to act accordingly.
        
    elif year == 2004:
        prod = pd.read_csv(r'G:\My Drive\k_data\CBS Households Expenditures Survey\\' + folder + '\\' + file_prod + '.txt', delim_whitespace = True, header = None)
        prod_pivoted = prod.pivot(index = 0, columns = 1, values = 2)
        mbs = pd.read_csv(r'G:\My Drive\k_data\CBS Households Expenditures Survey\\' + folder + '\\' + file_mbs + '.csv')
        
        # Renaming some of the columns to match the others.
        mbs.rename(columns = {'HOUSEHOLDS_NUMBER' : 'misparmb',
                              'N_PERSONS' : 'nefashot',
                              'WEIGHT' : 'weight',
                              'TYPE_OF_LOCALITY' : 'zurat_y'}, inplace = True)
        
        prod_pivoted.index = prod_pivoted.index.rename('misparmb')
 

        # Unfortunatly 2007 is also in a different format in the prod file, so yet again we load the file differently.        

    elif year == 2007:
        mbs = pd.read_csv(base_add + '\\' + folder + '\\' + file_mbs + '.csv')
        prod = pd.read_csv(base_add + '\\' + folder + '\\' + file_prod + '.txt', delim_whitespace = True, header = None)
        
        # Pivoting the long-formated prod file and renaming its index.
        
        prod_pivoted = prod.pivot(index = 1, columns = 2, values = 3)
        prod_pivoted.index = prod_pivoted.index.rename('misparmb')
        
        # Renaming the weights column for easier use later.
        
        mbs.columns = mbs.columns.str.lower()
        mbs.rename(columns = {'mishkal' : 'weight'}, inplace = True)


    # Finally, loading the files from 2008 and onwards.
    
    else:
        mbs = pd.read_csv(base_add + '\\' + folder + '\\' + file_mbs + '.csv')
        prod = pd.read_csv(base_add + '\\' + folder + '\\' + file_prod + '.csv')
        
        # Pivoting the long-formated pivot DataFrame.
        
        prod.columns = prod.columns.str.lower()
        prod_pivoted = prod.pivot(index = 'misparmb', columns = 'prodcode', values = 'schum')
        
        mbs.columns = mbs.columns.str.lower()
        
        # Renaming 2012-2014 weights column.
        
        mbs.rename(columns = {'mishkal' : 'weight'}, inplace = True)
    
    # Setting the household number as the index for easier joining columns fro the prod DataFrame.
    
    mbs.set_index('misparmb', inplace = True)    
    prod_pivoted.columns = prod_pivoted.columns.astype(str)
    
    # Joining some columns from the prod DataFrame to the household DataFrame: total gross income and income from allowances.
    
    mbs['i1_new'] = prod_pivoted['1'].fillna(0)
    mbs['i141_new'] = prod_pivoted['141'].fillna(0)
    mbs['i142_new'] = prod_pivoted['142'].fillna(0)    
    
    # Calculating gni per capita by substracting the allowances income from the total income, and then dividing the results by the households number of persons.
    
    mbs['gni'] = mbs['i1_new'] - mbs['i141_new'] - mbs['i142_new']
    mbs['gni_per_capita'] = mbs['gni'] / mbs['nefashot']
    
    # Creating a smal dictionary with the city filters for the households DataFrame.
    
    mask = dict(jeruz = mbs['zurat_y'] == 1,
                tlv = mbs['zurat_y'] == 2,
                israel = mbs['i1_new'] == mbs['i1_new'])
    
    # Calculating the average gni per capita, total residents, total gni and average gni per household.
    
    for city in mask.keys():
        results.loc[year, 'gni_' + city] = np.average(mbs.loc[mask[city], 'gni_per_capita'], weights = mbs.loc[mask[city], 'weight']) / cpi.loc[year, 'average'] * 100 * 12
        results.loc[year, 'capitas_' + city] = np.sum(mbs.loc[mask[city], 'nefashot'] * mbs.loc[mask[city], 'weight'])
        results.loc[year, 'total_gni' + city] = results.loc[year, 'gni_' + city] * results.loc[year, 'capitas_' + city]
        results.loc[year, 'gni_per_mb_' + city] = np.average(mbs.loc[mask[city], 'gni'], weights = mbs.loc[mask[city], 'weight']) / cpi.loc[year, 'average'] * 100 * 12

# Creating the plot.
    
anno = invert('מקור: סקרי הוצאות 9102-7991, הלמ"ס')

plt.figure(figsize = (10,5), dpi = 500)

plt.plot(results.index, results['gni_tlv'], label = invert('תל אביב'))
plt.plot(results.index, results['gni_israel'], label = invert('כלל המשק'))
plt.plot(results.index, results['gni_jeruz'], label = invert('ירושלים'))

plt.axvline(2011.5, linestyle = ':', label = invert('שינוי מתודולוגי'), color = 'black')
plt.axvline(2018.5, linestyle = '--', label = invert('שינוי מתודולוגי נוסף'), color = 'black')

plt.xticks(np.arange(period['start'] , period['end'] + 1), labels = list(map(str, np.arange(period['start'] , period['end'] + 1))), rotation = 90)
plt.yticks(plt.yticks()[0], ['{:,.0f}'.format(x) + '₪' for x in plt.yticks()[0]])
plt.title(invert('תל"ג לנפש בישראל ובערים שונות בשנים 9102-7991'), fontsize = 15)
plt.xlabel(invert('שנה'), fontsize = 12)
plt.ylabel(invert('ש"ח לנפש לשנה, מנוכה במדד המחירים לצרכן, בסיס 8102'), fontsize = 12)
plt.annotate(anno, (1995, 2000), annotation_clip = False)
plt.annotate('@tom_sadeh', (2016, 2000), annotation_clip = False)

plt.annotate('+' + str(np.round((results.loc[2019, 'gni_jeruz'] / results.loc[1997, 'gni_jeruz']-1) * 100, 1)) + '%', (2019, results.loc[2019, 'gni_jeruz'] - 300), color = 'tab:green')
plt.annotate('+' + str(np.round((results.loc[2019, 'gni_tlv'] / results.loc[1997, 'gni_tlv']-1) * 100, 1)) + '%', (2019, results.loc[2019, 'gni_tlv'] - 300), color = 'tab:blue')
plt.annotate('+' + str(np.round((results.loc[2019, 'gni_israel'] / results.loc[1997, 'gni_israel']-1) * 100, 1)) + '%', (2019, results.loc[2019, 'gni_israel'] - 300), color = 'tab:orange')

plt.legend(loc = 'upper left')
plt.grid(alpha = 0.5)
