from numpy import concatenate, zeros
import torch
from torch import nn
import numpy as np
from torch.autograd import Variable

LR = 0.02           # learning rate
class LSNN(nn.Module):
    def __init__(self, imput_siz, nodes_number, hidden_size = 10, hidden_layer = 1):
        super(LSNN, self).__init__()
        self.lstm = nn.LSTM(
            input_size = imput_siz,
            hidden_size = hidden_size,
            num_layers = hidden_layer,
            batch_first = True,
        )
        self.nodes_number = nodes_number
        self.hidden = (torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)),torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)))
        self.out = nn.Linear(hidden_size, 1)

    def forward(self,x):
        # x (batch, time_step, input_size)
        # h_state (n_layers, batch, hidden_size)
        # r_out (batch, time_step, output_size)
        r_out,self.hidden= self.lstm(x,self.hidden)
        self.hidden=(Variable(self.hidden[0]),Variable(self.hidden[1]))
        outs = []
        for time_step in range(self.nodes_number):
            outs.append(self.out(r_out[:, time_step, :]))
        return torch.stack(outs, dim=1)


def LSTM_train(pod_coeffs_all, history_level = 1):
    torch.manual_seed(1)    # reproducible
    print("LSTM_training")
    nodes_number = pod_coeffs_all.shape[1]
    time_steps = pod_coeffs_all.shape[0]
    lstmNN = LSNN(history_level, nodes_number, hidden_size = nodes_number)
    optimizer = torch.optim.Adam(lstmNN.parameters(), lr=LR)
    loss_func = nn.MSELoss()
    #prediction_list = []
    loss = []
    for step in range(time_steps - history_level):
        x_np = []
        y_np = []
        x = []
        y = []
        for i in range(history_level):
            x_np.append(pod_coeffs_all[step + i])
        x = Variable(torch.from_numpy(np.array(x_np)).float())
        y_np.append(pod_coeffs_all[step + history_level])
        y = Variable(torch.from_numpy(np.array(y_np)).float())
        x = x.permute(-1,0).view(1, nodes_number, history_level)
        y = y.permute(-1,0).view(1, nodes_number, 1)
        prediction = lstmNN(x)
        #prediction_list.append(prediction.data.view(nodes_number).numpy())
        loss = loss_func(prediction, y)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return lstmNN

def LSTM_predict(pod_coeffs_all, history_level, lstmNN, time_steps):
    print("LSTM_predicting")
    nodes_number = pod_coeffs_all.shape[1]
    prediction_list = []
    for step in range(time_steps - history_level):
        x_np = []
        y_np = []
        x = []
        y = []
        if(step < history_level):
            for i in range(history_level - step):
                x_np.append(pod_coeffs_all[step + i])
            for i in range(step + 1):
                if i > 0:
                    x_np.append(prediction_list[-i])
        elif(step >= history_level):
            for i in range(history_level):
                x_np.append(prediction_list[-history_level+i])
        x = Variable(torch.from_numpy(np.array(x_np)).float())
        x = x.permute(-1,0).view(1, nodes_number, history_level)
        with torch.no_grad():
            prediction = lstmNN(x)
        prediction_list.append(prediction.data.view(nodes_number).numpy())
        #print("length", len(prediction_list))
    return prediction_list
