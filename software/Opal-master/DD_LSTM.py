from numpy import concatenate, zeros
from scipy.linalg import toeplitz
import torch
from torch import nn
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from torch.autograd import Variable
import f2py


class LSNN(nn.Module):
    def __init__(self, local_nodes_number, neighbor_nodes_number, local_pos, hidden_size = 50, hidden_layer = 1):
        super(LSNN, self).__init__()
        self.lstm = nn.LSTM(  
            input_size = 1,
            hidden_size = hidden_size,    
            num_layers = hidden_layer,      
            batch_first = True,
        )
        self.local_nodes_number = local_nodes_number
        self.neighbor_nodes_number = neighbor_nodes_number
        self.local_pos = local_pos
        self.hidden = (torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)),torch.autograd.Variable(torch.zeros(hidden_layer, 1, hidden_size)))
        self.out = nn.Linear(hidden_size, 1)

    def forward(self,x):
        # x (batch, time_step, input_size)
        # h_state (n_layers, batch, hidden_size)
        # r_out (batch, time_step, output_size)
        r_out,self.hidden= self.lstm(x,self.hidden)
        self.hidden=(Variable(self.hidden[0]),Variable(self.hidden[1]))
        outs = []
        
        local_start = 0
        find_local = 0                                              # will turn to 1 if find local domain
        local_and_neighbor_nodes_number = 0
        for pos in range(len(self.neighbor_nodes_number)+1):         # number of neighbor and itself
            if(pos!=self.local_pos and find_local==0):              # have not find
                local_start += self.neighbor_nodes_number[pos]
            elif(pos==self.local_pos):                              # only happen once, find local domain
                find_local = 1
                local_start = local_start
                local_end = local_start + self.local_nodes_number
                local_and_neighbor_nodes_number = local_end
            else:                                                   # after find
                local_and_neighbor_nodes_number += self.neighbor_nodes_number[pos - find_local]

        
        for time_step in range(local_and_neighbor_nodes_number):
            if(time_step>=local_start and time_step<local_end):
                outs.append(self.out(r_out[:, time_step, :]))
        return torch.stack(outs, dim=1)

def DD_LSTM routine(
    LR = 0.02           # learning rate

    nx = 128
    ny = 1
    nsplt = 3
    whichd, ncola, cola, fina = my.DD(nx, ny, nsplt)

    print(whichd)


    NSPLIT = pow(2, nsplt)

    cell = nx/NSPLIT

    nodes = np.arange(nx)

    list_subs_all = []

    for sub in range(NSPLIT):                             # loop over all subdomains

        list_subs = [sub]                                 # include subdomain sub in list

        for node in nodes:                               # loop over all  nodes

            if(whichd[node] == sub):                    # look for colns that belong to another subdomain

                for count in range(int(fina[node]), int(fina[node + 1])):         # look at non-zeros in row nod

                    col = cola[count]               # this is the coln of the fem graph
                #print(col)
                    sub2 = whichd[col]

                    if sub2 != sub:                    # add to the list of subdomains connected

                        list_subs.append(sub2)

        list_subs_all.append(list_subs)

    print "list_subs_all", list_subs_all


    reorder = my.reorde(list_subs_all, NSPLIT)



    torch.manual_seed(1)    # reproducible
    mat.use("TkAgg")


    local_nodes_number = []
    neighbor_nodes_number = []
    for i in range(NSPLIT):
        local_nodes_number.append(whichd.tolist().count(i))
        one_neighbor_nodes_number = []
        for j in list_subs_all[i]:                                 # self and neighbor
            if j != i:
	        #print("j = ", j, "   i = ", i)
	        one_neighbor_nodes_number.append(whichd.tolist().count(j))
	        #print("one_neighbor_nodes_number", one_neighbor_nodes_number)
        neighbor_nodes_number.append(one_neighbor_nodes_number)
        
    print("local node number", local_nodes_number)
    print("neighbor_nodes_number", neighbor_nodes_number)


    local_pos = np.zeros(NSPLIT).tolist()

    total_domain_number = NSPLIT
    lstmNN = []
    optimizer = []

    for i in range(total_domain_number):
        lstmNN.append(LSNN(local_nodes_number[i], neighbor_nodes_number[i], local_pos[i]))
        optimizer.append(torch.optim.Adam(lstmNN[i].parameters(), lr=LR))

    
    loss_func = nn.MSELoss()

    prediction = np.zeros(NSPLIT).tolist()
    loss = np.zeros(NSPLIT).tolist()

    steps = np.linspace(0, nx, nx, dtype=np.float32)

    for step in range(126):
        x_np = []
        y_np = []
        x = []
        y = []
        if step == 0:
	        for sub in range(total_domain_number):
                x_np.append(toeplitz(concatenate([[1.], zeros(127)]),concatenate([[1.,1.,1.], zeros(125)]))[step : step + 1, cell*sub:cell*(sub+1)])
                y_np.append(toeplitz(concatenate([[1.], zeros(127)]),concatenate([[1.,1.,1.], zeros(125)]))[step + 1: step + 2, cell*sub:cell*(sub+1)])
	        for sub in range(total_domain_number):
	            sub_list = []
	            for local_sub in list_subs_all[sub]:
	                sub_list.append(x_np[local_sub])
	            x.append(Variable(torch.from_numpy(np.array(sub_list)).float()))
	            y.append(Variable(torch.from_numpy(np.array(y_np[sub])).float()))
        else:
	        for sub in range(total_domain_number):
                y_np.append(toeplitz(concatenate([[1.], zeros(127)]),concatenate([[1.,1.,1.], zeros(125)]))[step + 1: step + 2, cell*sub:cell*(sub+1)])

	    for sub in range(total_domain_number):
	        sub_list = [prediction[sub].data.view(local_nodes_number[sub]).numpy().tolist(), ]           # itself's nodes number
	        for local_sub in list_subs_all[sub][1:]:
	            sub_list.append(prediction[local_sub].data.view(cell).numpy().tolist())   # it neighbor's node number
	        x.append(Variable(torch.from_numpy(np.array(sub_list)).float()))
	        y.append(Variable(torch.from_numpy(np.array(y_np[sub])).float()))

    
        #print("y[0].size()", y[1].size())

        for i in range(total_domain_number):
	    x[i] = x[i].view(1,-1,1)
	    y[i] = y[i].view(1,-1,1)
        #print("x_type", type(x), x.size())

        for i in range(total_domain_number):
            prediction[i] = lstmNN[i](x[i])
            loss[i] = loss_func(prediction[i], y[i])


        for i in range(total_domain_number):
            optimizer[i].zero_grad()
        
        for i in range(total_domain_number):  # backpropagation, compute gradients
            loss[i].backward()

    
        for i in range(total_domain_number):
            optimizer[i].step()
        #print(loss_func(prediction, x1[:,:33,:]).detach().numpy())
        for i in range(20):
        #while(loss_func(prediction1, x1[:,:33,:]) > 0.01 or loss_func(prediction2, x2[:,33: 66,:])>0.01 or loss_func(prediction3, x3[:,66: ,:])>0.01):
	    x_tem = []
	    y_tem = []
	    for sub in range(total_domain_number):
	        #print(x[sub][:,:local_nodes_number[sub],:].data.view(local_nodes_number[sub]).numpy().tolist())
	        sub_list = [x[sub][:,:local_nodes_number[sub],:].data.view(local_nodes_number[sub]).numpy().tolist(), ]           # itself's nodes number
	        for local_sub in list_subs_all[sub][1:]:
		    #print("type", local_sub_node_num,sub)
	            sub_list.append(prediction[local_sub].data.view(cell).numpy().tolist())   # it neighbor's node number
	        x_tem.append(Variable(torch.from_numpy(np.array(sub_list)).float()))
	        y_tem.append(Variable(torch.from_numpy(np.array(y_np[sub])).float()))

    	for i in range(total_domain_number):
	        x[i] = x_tem[i].view(1,-1,1)
	        y[i] = y_tem[i].view(1,-1,1)
        
        
        for i in range(total_domain_number):
            prediction[i] = lstmNN[i](x[i])
            loss[i] = loss_func(prediction[i], y[i])
        
        for i in range(total_domain_number):
            optimizer[i].zero_grad()
        
        for i in range(total_domain_number):  # backpropagation, compute gradients
            loss[i].backward()
            
        for i in range(total_domain_number):
            optimizer[i].step()
    
        y_show = []
        prediction_show = []
        for i in range(total_domain_number):
	        y_show.append(y[i].data.numpy().flatten())
	        prediction_show.append(prediction[i].data.numpy().flatten())

    
        plt.figsize=(20, 10)
        plt.ion()
        plt.title(step,fontsize=24)
    
        plt.plot(steps, np.array(y_show).reshape(nx, ny), 'r-')
        plt.plot(steps, np.array(prediction_show).reshape(nx, ny), 'b-')
        plt.legend()
        plt.draw()
        plt.pause(0.01)
        plt.clf()

    plt.show()
