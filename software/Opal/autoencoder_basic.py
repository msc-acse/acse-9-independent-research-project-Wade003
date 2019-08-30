#!/usr/bin/env python

#    Copyright (C) 2017 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

#Authors: T. Phillips, C.E. Heaney and C.C. Pain 

import os
os.environ['KERAS_BACKEND'] = 'tensorflow'
import keras
from keras.layers import Input, Dense, Lambda
from keras.models import Model
from keras.callbacks import TensorBoard
from keras import regularizers



#print "keras:", os.path.abspath(keras.__file__)

# missing from activations.py file, so have changed elu to relu
#def elu(x, alpha=1.):
#    return K.elu(x, alpha)

# Basic Model 
#each row of a matrix contains each snapshot used in the model
#with each column containing a variable

# inputs
# number of snapshots
# i = 1000 # now m
# full number of variables
# n = 500  # still n
# Data
#'''U_train = np.zeros((i,n))''' now U_train((m,n))

# reduced number
#nr = 100

def autoencoder(m,n,U_train,nr,nirom_options):

   my_batch_size = nirom_options.compression.batch_size 
   nepochs =  nirom_options.compression.nepochs

   print "my batch size", my_batch_size
   print "nepochs", nepochs

   U_test = U_train
   # input data
   input_img = Input(shape=(n,))
   #encoding layers # 382 320 255 192
   encoded = Dense(382, activation='elu')(input_img)
   encoded = Dense(320,activation='elu')(encoded)
   encoded = Dense(255,activation='elu')(encoded)
   encoded = Dense(192,activation='elu')(encoded)
   encoded = Dense(nr,activation='elu')(encoded)

   #decoding layers
   decoded = Dense(192,activation='elu')(encoded)
   decoded = Dense(255,activation='elu')(decoded)
   decoded = Dense(320,activation='elu')(decoded)
   decoded = Dense(382,activation='elu')(decoded)
   decoded = Dense(n,activation='elu')(decoded)
   #putting the model together (input,last layer)
   autoencoder=Model(input_img,decoded)
   #determine optimizer and loss functiosn
   autoencoder.compile(optimizer='Nadam',loss='mean_squared_error')

   #train the model. (input data, output data (identical for autoencoder),
   # iterations, batch size, shuffle data randomly, test data)
   '''history = autoencoder.fit(U_train,U_train,epochs = 40000,batch_size = t,
                #shuffle = True, validation_data=(U_test, U_test),)'''
   history = autoencoder.fit(U_train,U_train,epochs = nepochs,batch_size = my_batch_size,
                shuffle = True, validation_data=(U_test, U_test),)

   # retrieve encoder section
   input_img = Input(shape=(n,)) #382 320 255 192
   encoded = Dense(382,activation='elu',
                weights=autoencoder.layers[1].get_weights())(input_img)
   encoded = Dense(320,activation='elu',
                weights=autoencoder.layers[2].get_weights())(encoded)
   encoded = Dense(255,activation='elu',
                weights=autoencoder.layers[3].get_weights())(encoded)
   encoded = Dense(192,activation='elu',
                weights=autoencoder.layers[4].get_weights())(encoded)
   encoded = Dense(nr,activation='elu',
                weights=autoencoder.layers[5].get_weights())(encoded)
   encoder = Model(input_img, encoded)
   # retrieve decoder section
   encoded_input = Input(shape=(nr,))
   decoded_out = Dense(192,activation='elu',
                    weights=autoencoder.layers[6].get_weights())(encoded_input)
   decoded_out = Dense(255,activation='elu',
                    weights=autoencoder.layers[7].get_weights())(decoded_out)
   decoded_out = Dense(320,activation='elu',
                    weights=autoencoder.layers[8].get_weights())(decoded_out)
   decoded_out = Dense(382,activation='elu',
                    weights=autoencoder.layers[9].get_weights())(decoded_out)
   decoded_out = Dense(n,activation='elu',
                    weights=autoencoder.layers[10].get_weights())(decoded_out)
   decoder = Model(encoded_input, decoded_out)

   # Used to predict results
   U_predict = autoencoder.predict(U_train) # encoder.predict, decoder.predict

   print "U_train.shape, U_predict.shape", U_train.shape, U_predict.shape
   # Used to save model for future use
   encoder.save('encoder.h5')
   decoder.save('decoder.h5')

   U_encoded = encoder.predict(U_train)
   print "shape U_encoded", U_encoded.shape 
   # Used to load model in another script. Encoder and Decoder can be retreved from this model. 
   '''autoencoder = load_model('autoencoder.h5')'''
   autoencoder.summary()

   return U_predict, U_encoded


