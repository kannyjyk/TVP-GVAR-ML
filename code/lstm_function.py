# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 11:08:13 2021

@author: zhang
"""

import numpy as np
from keras.models import Sequential
from keras.layers import Dense, LSTM, Dropout
from keras.models import Model
from keras.layers import Input, concatenate
from keras.layers.recurrent import GRU
import warnings
warnings.filterwarnings('ignore')

def create_dataset(dataset, look_back=1):
    dataX, dataY = [], []
    for i in range(len(dataset)-look_back):
        dataX.append(dataset[i:(i+look_back)])
        dataY.append(dataset[i + look_back])
    return np.array(dataX), np.array(dataY)

def lstm_gvar(bt_train, h_pred = 12):
    np.random.seed(2021)
    trainX, trainY = create_dataset(bt_train)
    trainX = np.reshape(trainX, (trainX.shape[0], 1, trainX.shape[2]))

    batch_size = 1
    epochs = 20
    model = Sequential()
    model.add(LSTM(32, 
                    batch_input_shape=(batch_size, 1, trainX.shape[2]),
                    stateful=True))
    model.add(Dense(trainY.shape[1]))
    model.compile(loss='mean_squared_error', optimizer="adam", metrics=["accuracy"])    
    model.fit(trainX, trainY, epochs=epochs, batch_size=batch_size, verbose=0 , shuffle=False)
    
    x = np.vstack((trainX[-1][1:], (trainY[-1])))
    preds = []
    pred_num = h_pred
    for i in np.arange(pred_num):
        pred = model.predict(x.reshape((x.shape[0], 1, x.shape[1])), batch_size = batch_size)
        preds.append(pred.squeeze())
        x = np.vstack((x[1:],pred))
    
    preds = np.array(preds)
    return preds

def gru_gvar(bt_train, h_pred = 12):
    np.random.seed(2021)
    trainX, trainY = create_dataset(bt_train)
    trainX = np.reshape(trainX, (trainX.shape[0], 1, trainX.shape[2]))

    batch_size = 1
    epochs = 20
    model = Sequential()
    model.add(GRU(32, 
                  batch_input_shape=(batch_size, 1, trainX.shape[2]),
                  stateful=True))
    model.add(Dense(trainY.shape[1]))
    model.compile(loss='mean_squared_error', optimizer="adam", metrics=["accuracy"])    
    model.fit(trainX, trainY, epochs=epochs, batch_size=batch_size, verbose=0 , shuffle=False)
    
    x = np.vstack((trainX[-1][1:], (trainY[-1])))
    preds = []
    pred_num = h_pred
    for i in np.arange(pred_num):
        pred = model.predict(x.reshape((x.shape[0], 1, x.shape[1])), batch_size = batch_size)
        preds.append(pred.squeeze())
        x = np.vstack((x[1:],pred))
    
    preds = np.array(preds)
    return preds

def lstm_gru_gvar(bt_train, h_pred = 12):
    np.random.seed(2021)
    trainX, trainY = create_dataset(bt_train)
    trainX = np.reshape(trainX, (trainX.shape[0], 1, trainX.shape[2]))

    batch_size = 1
    epochs = 20
    
    inputs = Input(batch_shape = (batch_size, 1, trainX.shape[2]))
    x1 = GRU(32, batch_input_shape=(batch_size, 1, trainX.shape[2]), stateful=True)(inputs)
    x2 = LSTM(32, batch_input_shape=(batch_size, 1, trainX.shape[2]), stateful=True)(inputs)
    x3 = concatenate([x1, x2], axis = 1)
    
    x4 = Dense(trainY.shape[1])(x3)
    model = Model(inputs, x4)

    model.compile(loss='mean_squared_error', optimizer="adam", metrics=["accuracy"])    
    model.fit(trainX, trainY, epochs=epochs, batch_size=batch_size, verbose=0 , shuffle=False)
    
    x = np.vstack((trainX[-1][1:], (trainY[-1])))
    preds = []
    pred_num = h_pred
    for i in np.arange(pred_num):
        pred = model.predict(x.reshape((x.shape[0], 1, x.shape[1])), batch_size = batch_size)
        preds.append(pred.squeeze())
        x = np.vstack((x[1:],pred))
    
    preds = np.array(preds)
    return preds
