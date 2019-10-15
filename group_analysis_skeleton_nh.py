#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scene-cat problem set for PSY 1210 - Fall 2018

@author: Michael Mack
"""

#%% import block 
import numpy as np
import scipy as sp
import scipy.stats
import os
import shutil


#%%
# Via sp.savetxt, copy files from testing room folders to raw data, rename files to include
# testing room letter in the filename
#also include my initials 

# to create an empty array in NumPy to store my matrix
data=np.empty((0,5))

testingrooms = ['A','B','C']
for room in testingrooms:
   file_x = "testingroom"+room+"/experiment_data.csv"
   tmp = sp.loadtxt(file_x, delimiter=',')
   data=np.vstack([data,tmp])
   sp.savetxt('rawdata/'+'experimentNH'+room+'.csv',tmp)
  



#%%
# read in all the data files in rawdata directory using a for loop 
#also read in data via sp.loadtxt and np.vstack
# columns: subject, stimulus, pairing, accuracy, median RT

# 5 columns of Subject, Stimuli, Pairing, Accuracy, Median RT
data=np.empty((0,5))

testingrooms = ['A','B','C']
for room in testingrooms:
   file_x = "testingroom"+room+"/experiment_data.csv"
   tmp = sp.loadtxt(file_x, delimiter=',')
   data=np.vstack([data,tmp])
   sp.savetxt('rawdata/'+'experimentNH'+room+'.csv',tmp)


#%%
# calculate overall average accuracy and average median RT

# set acc as new object of list
# Specify element range of 92 as there are 92 participants
# Set x as new object of data, element, and 3rd row of accuracy
# Append x to acc list 
acc = []
for element in range(92):
    x = data[element][3]
    acc.append(x)
    

# Calculate the mean of accuracy via numpy's mean function for acc_avg
acc_avg = np.mean(acc)

#acc_avg = np.mean(acc)*100...   # 91.48%


# Apply similar steps for above accuracy calculations to mrt 
mrt = []
for element in range(92):
    y = data[element][4]
    mrt.append(y)
    
mrt_avg = np.mean(mrt)

#mrt_avg = ...   # 477.3ms
#%%
# calculate averages (accuracy & RT) split by stimulus using a for loop and an 
# if statement. (i.e., loop through the data to make a sum for each condition, 
# then divide by the number of data points going into the sum)

# Create 4 new lists
acc_words = []
mrt_words = []
acc_faces = []
mrt_faces = []

# Set z to equal to 1 which equates to stimulus value for words 
# Then append accuracy from 3rd column and median RT from 4th column to correct lists
for element in range(92):
    z = data[element][1]
    if z==1:
        acc_words.append(data[element][3])
        mrt_words.append(data[element][4])
        
    else:
        acc_faces.append(data[element][3])
        mrt_faces.append(data[element][4])

        
# Calculate mean accuracy and median RTs for word and face stimuli separately         
acc_words_avg = np.mean(acc_words)
mrt_words_avg = np.mean(mrt_words)
acc_faces_avg = np.mean(acc_faces)
mrt_faces_avg = np.mean(mrt_faces)

# words: 88.6%, 489.4ms   faces: 94.4%, 465.3ms


#%%
# calculate averages (accuracy & RT) split by congruency using indexing, 
# slicing, and numpy's mean function 
# wp - white/pleasant, bp - black/pleasant
# (hint: only one line of code is needed per average)
# Similar to above steps for Stimuli, set v to equal to 1 which equates to pairing value of white_pleasant
# Then append accuracy from 3rd column and median RT from 4th column to correct lists

acc_wp = []
mrt_wp = []
acc_bp = []
mrt_bp = []

for element in range(92):
    v = data[element][2]
    if v==1:
        acc_wp.append(data[element][3])
        mrt_wp.append(data[element][4])
        
    else:
        acc_bp.append(data[element][3])
        mrt_bp.append(data[element][4])
        
acc_wp_avg = np.mean(acc_wp)
mrt_wp_avg = np.mean(mrt_wp)
acc_bp_avg = np.mean(acc_bp)
mrt_bp_avg = np.mean(mrt_bp)

acc_wp = ...  # 94.0%
acc_bp = ...  # 88.9%
mrt_wp = ...  # 469.6ms
mrt_bp = ...  # 485.1ms


#%% 
# calculate average median RT for each of the four conditions
# use for loops, indexing/slicing, or both!
# (hint: might be easier to slice data into separate words and faces datasets)


# Create 4 new lists for the four conditions
mrt_words_wp = []
mrt_words_bp = []
mrt_faces_wp = []
mrt_faces_bp = []

# For stimuli in 1st column, set w to equal to 1 which equates to stimulus value for words 
# For pairing in 2nd column, set v to equal to 1 which equates to pairing value of white_pleasant
for element in range(92):
    w = data[element][1]
    v = data[element][2]
    if w == 1 and v == 1:
        mrt_words_wp.append(data[element][4])
    if w == 1 and v ==2:
        mrt_words_bp.append(data[element][4])
    if w == 2 and v ==1:
        mrt_faces_wp.append(data[element][4])   
    if w == 2 and v ==2:
        mrt_faces_bp.append(data[element][4])

mrt_words_wp_avg = np.mean(mrt_words_wp)
mrt_words_bp_avg = np.mean(mrt_words_bp)
mrt_faces_wp_avg = np.mean(mrt_faces_wp)
mrt_faces_bp_avg = np.mean(mrt_faces_bp)

# words - white/pleasant: 478.4ms
# words - black/pleasant: 500.3ms
# faces - white/pleasant: 460.8ms
# faces - black/pleasant: 469.9ms


#%%        
# compare pairing conditions' effect on RT within stimulus using scipy's 
# paired-sample t-test: scipy.stats.ttest_rel()
#
import scipy.stats

# Extract t-value from scipy result via function of .statistic
# Extract p-value from scipy result via function of .pvalue
ttest_word = scipy.stats.ttest_rel(mrt_words_wp, mrt_words_bp).statistic
pvalue_word = scipy.stats.ttest_rel(mrt_words_wp, mrt_words_bp).pvalue
ttest_face = scipy.stats.ttest_rel(mrt_faces_wp, mrt_faces_bp).statistic
pvalue_face = scipy.stats.ttest_rel(mrt_faces_wp, mrt_faces_bp).pvalue

# words: t=-5.36, p=2.19e-5
# faces: t=-2.84, p=0.0096


#%%
# print out averages and t-test results
# (hint: use the ''.format() method to create formatted strings)

#Overall averages
print('\nOVERALL_Percent_Accuracy_and_Median_Reaction_Time: {:.2f}%, {:.1f} ms'.format(100*acc_avg,mrt_avg))

#Split by stimulus
print('\nWORD_STIMULI_Percent_Accuracy_and_Median_Reaction_Time: {:.2f}%, {:.1f} ms'.format(100*acc_words_avg, mrt_words_avg))
print('\nFACE_STIMULI_Percent_Accuracy_and_Median_Reaction_Time: {:.2f}%, {:.1f} ms'.format(100*acc_faces_avg, mrt_faces_avg))

#Split by congruency
print('\nCONGRUENCY(WHITE/PLEASANT)_Percent_Accuracy_and_Median_Reaction_Time: {:.2f}%, {:.1f} ms'.format(100*acc_wp_avg, mrt_wp_avg))
print('\nCONGRUENCY(BLACK/PLEASANT)_Percent_Accuracy_and_Median_Reaction_Time: {:.2f}%, {:.1f} ms'.format(100*acc_bp_avg, mrt_bp_avg))

#Average median RT for each of the four conditions
print('\nWORD(WHITE/PLEASANT)_Median_Reaction_Time: {:.1f} ms'.format(mrt_words_wp_avg))
print('\nWORD(BLACK/PLEASANT)_Median_Reaction_Time: {:.1f} ms'.format(mrt_words_bp_avg))
print('\nFACE(WHITE/PLEASANT)_Median_Reaction_Time: {:.1f} ms'.format(mrt_faces_wp_avg))
print('\nFACE(BLACK/PLEASANT)_Median_Reaction_Time: {:.1f} ms'.format(mrt_faces_bp_avg))

#ttest and pvalue for comparing median RT for the two pairings for each stimulus 
print('\nt_Value_WORD: {:.2f}'.format(ttest_word))
print('\np_Value_WORD: {:.5f}'.format(pvalue_word))
print('\nt_Value_FACE: {:.2f}'.format(ttest_face))
print('\np_Value_FACE: {:.5f}'.format(pvalue_face))


