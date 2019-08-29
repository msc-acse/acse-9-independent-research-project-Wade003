from numpy import concatenate, zeros
from scipy.linalg import toeplitz
import matplotlib.pyplot as plt

import torch
import torch.nn as nn

from numpy import concatenate, zeros
from scipy.linalg import toeplitz
import torch
from torch import nn
import numpy as np
import matplotlib as mat
mat.use("TkAgg")
import matplotlib.pyplot as plt
import time
from torch.autograd import Variable
import cv2

torch.manual_seed(1)    # reproducible

hidden_siz = 50
hidden_lay = 1

LR = 0.02           # learning rate
class LSNN(nn.Module):
    def __init__(self):
        super(LSNN, self).__init__()
        self.lstm = nn.LSTM(  
            input_size=5,
            hidden_size=hidden_siz,    
            num_layers=hidden_lay,      
            batch_first=True,

        )
        self.hidden = (torch.autograd.Variable(torch.zeros(hidden_lay, 1, hidden_siz)),torch.autograd.Variable(torch.zeros(hidden_lay, 1, hidden_siz)))
        self.out = nn.Linear(hidden_siz, 1)
        
    def forward(self,x):
        # x (batch, time_step, input_size)
        # h_state (n_layers, batch, hidden_size)
        # r_out (batch, time_step, output_size)
        r_out,self.hidden= self.lstm(x,self.hidden)

        self.hidden=(Variable(self.hidden[0]),Variable(self.hidden[1]))
        outs = []
        for time_step in range(r_out.size(1)):
            outs.append(self.out(r_out[:, time_step, :]))
        return torch.stack(outs, dim=1)


lstmNN = LSNN()
optimizer = torch.optim.Adam(lstmNN.parameters(), lr=LR)  # optimize all rnn parameters
loss_func = nn.MSELoss()

loss_list = []
prediction_list = []
for step in range(80-6):
    steps = np.linspace(0, 100, 100, dtype=np.float32)
    
    if step == 0:
        x_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step: 5, :]
        y_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[5:6, :] 
    else:
        x_np = concatenate([
                            toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step: step+4, :],
                            prediction.view(1,100).data.numpy()])
        y_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step+5:step+6, :]
    
    #x_np = steps    # float32 for converting torch FloatTensor
    #y_np = steps
    x = Variable(torch.from_numpy(x_np).float())  # shape (batch, time_step, input_size)
    y = Variable(torch.from_numpy(y_np).float())
    x = x.permute(-1,0).view(1,100,5)
    y = y.view(1,100,1)
    prediction = lstmNN(x)
    #print(prediction.size())
    prediction_list.append(prediction.data.view(100).numpy())

    loss = loss_func(prediction, y)     # cross entropy loss
    loss_list.append(loss)
    
    #train_loss += loss*X.size(0)
    
    optimizer.zero_grad()               # clear gradients for this training step
    loss.backward()                     # backpropagation, compute gradients
    optimizer.step()
    # apply gradients
    plt.figsize=(20, 10)
    plt.ion()
    plt.title(step,fontsize=24)
    plt.plot(steps, y_np.flatten(), 'r-')
    plt.plot(steps, prediction.data.numpy().flatten(), 'b-')
    #plt.legend()
    plt.draw()
    plt.pause(0.1)
    plt.clf()


plt.show()

plt.plot(steps[:80-6], loss_list, label = 'Loss')
plt.legend()
plt.show()

torch.manual_seed(1)    # reproducible

loss_func = nn.MSELoss()

loss_list = []
prediction_list = []
for step in range(100-6):
    steps = np.linspace(0, 100, 100, dtype=np.float32)
    
    if step <6:
        if step == 0:
            x_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step: 5, :]
            y_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[5:6, :] 
        else:
            
            x_np = concatenate([toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step: 5, :],
                                np.array(prediction_list[:step])])
            y_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step+5:step+6, :]
        
    else:
        x_np = np.array(prediction_list[-5:])
        y_np = toeplitz(concatenate([[1.], zeros(99)]),concatenate([[1.,1.,1.], zeros(97)]))[step+5:step+6, :]   
    
    #x_np = steps    # float32 for converting torch FloatTensor
    #y_np = steps
    x = Variable(torch.from_numpy(x_np).float())  # shape (batch, time_step, input_size)
    y = Variable(torch.from_numpy(y_np).float())

    x = x.permute(-1,0).view(1,100,5)
    y = y.view(1,100,1)
    with torch.no_grad():
        prediction = lstmNN(x)
        #print("pre ",prediction.data.size())
        prediction_list.append(prediction.data.view(100).numpy())
        #print(prediction_list)
        loss = loss_func(prediction, y)     # cross entropy loss
        loss_list.append(loss)
    
    #train_loss += loss*X.size(0)
    
    # apply gradients
    plt.figsize=(20, 10)
    plt.ion()
    plt.title(step,fontsize=24)
    plt.plot(steps, y_np.flatten(), 'r-')
    plt.plot(steps, prediction.data.numpy().flatten(), 'b-')
    #plt.legend()
    plt.draw()
    plt.pause(0.1)
    plt.clf()

    plt.ioff()
    #plt.show()
plt.plot(steps, y_np.flatten(), 'r-')
plt.plot(steps, prediction.data.numpy().flatten(), 'b-')
plt.show()

plt.plot(steps[:94], loss_list, label = 'Error')
plt.legend()
plt.show()