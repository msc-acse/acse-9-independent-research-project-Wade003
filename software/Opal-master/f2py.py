import decomp_new
import numpy as np
import vtktools

def get_coordinates(vtu_filename):
    vtu_data = vtktools.vtu(vtu_filename)
    coordinates = vtu_data.GetLocations()
    return vtu_data,coordinates

def weight_from_stress_tensor(snapshots, nsnapshot, ndim, nonods):
    a=0.05
    b=0.0019183
    vel_snap = np.zeros((ndim, nsnapshot))
    vel_mean = np.zeros(ndim)
    reynolds =  np.zeros(ndim)
    wnod = np.zeros(nonods)
    for node in range(nonods):
        for i in range(ndim):
            vel_snap[i, :] = snapshots[:, i*nonods + node]
            vel_mean[i] = sum(vel_snap[i,:])/float(nsnapshot)
            reynolds[i] = sum((vel_snap[i,:] - vel_mean[i])**2 )/float(nsnapshot)
        max_reynolds = max(reynolds)
        wnod[node] = (1./b) * np.log( max_reynolds/a + 1.0)
    return wnod

def DD(nNodes, nsplt,fwd_options, nirom_options):
    #print(decomp.python_set_up_recbis.__doc__)

    #n = nx*ny
    #ncola, cola, fina = decomp_new.one_d_row_stor(nNodes,1,nNodes,(nNodes*3)-2)
    #ncola, cola, fina = decomp_new.one_d_row_stor(nx,ny,n,(n*5))

    ext = fwd_options.results_extension
    filename = fwd_options.path_to_results + '/' + fwd_options.results_filebase + '_0.' + ext
    vtu_data, coordinates = get_coordinates(filename)
    print("coordinates", type(coordinates), coordinates.shape)
    #wnod = decomp_new.weight_from_stress_tensor(nsplt, ndim, nonods)
    findm,colm,ncolm = decomp_new.fsmfp(coordinates, nNodes)
    fina = findm
    cola = colm
    ncola = ncolm
    print("fina", fina)
    split_levels = np.zeros((nsplt),dtype='int32')
    split_levels[:] = 2

    havwnod = 2
    #print("nirom_options.snapshots.values[0]", nirom_options.snapshots.values[0].transpose(1, 0).shape, nNodes)

    wnod = np.zeros(nNodes)
    #wnod = weight_from_stress_tensor(nirom_options.snapshots.values[0].transpose(1, 0), nirom_options.snapshots.values[0].shape[1], fwd_options.ndim, nNodes)

    # wnod = decomp_new.weight_from_stress_tensor(nirom_options.snapshots.values[0].transpose(1, 0), nirom_options.snapshots.values[0].shape[1], fwd_options.ndim, nNodes)
    #
    wnod[:] = 1.0

    havmat = 0
    a = np.zeros(1)
    exact = True
    iexact = 1
    ii=1
    na=0

    print "OK up to here"

    #witchd = decomp_new.python_set_up_recbis(splev,fina,cola, nsplt,ncola,n, havwnod,wnod,exact, havmat)

    print "split_levels", split_levels
    print "shape fina", fina.shape
    print "shape cola", cola.shape
    print "nsplit", nsplt
    print "ncola", ncola
    print "ii", ii
    print("nNodes", nNodes)
    #print decomp_new.python_set_up_recbis.__doc__
    # we have to finish the sub clal with the variables used in the decleration within this sub.
    print("about to print whichd")
    witchd = decomp_new.python_set_up_recbis(split_levels,fina,cola, wnod,a, havwnod,havmat,iexact, nsplt,ncola,nNodes,na)
    print(witchd - 1)
    print "ncola", ncola
    print "cola", cola[:ncola] - 1, len(cola)
    print "fina", fina - 1, len(fina)
    witchd =  witchd - 1
    #visit_list = []
    #for i in range(len(witchd)):
	#if witchd[i] not in visit_list:
	 #   visit_list.append(witchd[i])
    #item_number = []
    #for item in visit_list:
    #    item_number.append(witchd.tolist().count(item))
    #new_witchd = []
    #for i in range(len(visit_list)):
    #    for j in range(item_number[i]):
#	    new_witchd.append(len(visit_list) -1 - i)
	    #new_witchd.append(i)

    return np.array(witchd), ncola, cola[:ncola] - 1, fina - 1 #

def get_local_list(whichd, fina, cola, nsplt, node_number):
    NSPLIT = pow(2, nsplt)

    cell = node_number/NSPLIT

    nodes = np.arange(node_number)

    local_list = []

    for sub in range(NSPLIT):                             # loop over all subdomains
        list_subs = [sub]                                 # include subdomain sub in list

        for node in nodes:                               # loop over all  nodes

            if(whichd[node] == sub):                    # look for colns that belong to another subdomain

                for count in range(int(fina[node]), int(fina[node + 1])):         # look at non-zeros in row nod

                    col = cola[count]               # this is the coln of the fem graph
                    #print(col)
                    sub2 = whichd[col]

                    if sub2 != sub and sub2 not in list_subs:                    # add to the list of subdomains connected

                        list_subs.append(sub2)

        local_list.append(list_subs)
    return local_list



def reorde(local_list, num_sub):
    count_fina = 1
    fina = [count_fina]
    cola = []
    for sub_and_neighbor in local_list:
        count_fina = count_fina + len(sub_and_neighbor)
        fina.append(count_fina)
        for item in sub_and_neighbor:
            cola.append(item+1)
    print("SUB_fina", fina)
    print("SUB_cola", cola)
    print("num_sub", num_sub)
    order = decomp_new.reorde(fina,cola, num_sub,len(cola))
    return order - 1

def LSTM_train(dd_pod_coeffs_all, history_level, local_list):
    torch.manual_seed(1)    # reproducible
    print("DD_LSTM_training")
    lstmNN = []
    optimizer = []
    time_steps = dd_pod_coeffs_all.shape[1]
    nodes_number = dd_pod_coeffs_all.shape[2]
    for sub in range(len(local_list)):
        lstmNN.append(LSNN(history_level, nodes_number))
        optimizer.append(torch.optim.Adam(lstmNN[sub].parameters(), lr=LR))
    loss_func = nn.MSELoss()
    loss_1 = []
    n_sub = len(local_list)
    for step in range(time_steps - history_level):
        # prediction_list = []
        # for sub in range(n_sub):
        #     x_np = []
        #     y_np = []
        #     x = []
        #     y = []
        #     for i in range(history_level):
        #         null_list = []
        #         for local_input in local_list[sub]:
        #             null_list.append(dd_pod_coeffs_all[local_input][step + i])
        #         # null_list.append(dd_pod_coeffs_all[local_list[sub][0]][step + i])
        #         x_np.append(null_list)
        #     x = Variable(torch.from_numpy(np.array(x_np)).float())
        #     y_np.append(dd_pod_coeffs_all[sub][step + history_level])
        #     y = Variable(torch.from_numpy(np.array(y_np)).float())
        #     x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
        #     # x = x.permute(1,2,0).view(-1, nodes_number, history_level)
        #     y = y.permute(-1,0).view(-1, nodes_number, 1)
        #
        #     # x = x.permute(1,0).view(1, nodes_number*len(local_list[sub]), history_level)
        #     # y = y.permute(1,0).view(1, nodes_number, 1)
        #     prediction = lstmNN[sub](x)
        #     prediction_list.append(prediction.data.view(nodes_number).numpy())
        #     loss = loss_func(prediction, y)
        #     optimizer[sub].zero_grad()
        #     loss.backward()
        #     optimizer[sub].step()

        prediction_update = []
        if(step == 0):
            for sub in range(len(local_list)):
                prediction_update.append([])
                prediction_update[sub].append(dd_pod_coeffs_all[sub][step + history_level - 1])
        else:
            for sub in range(len(local_list)):
                prediction_update.append([])
                prediction_update[sub].append(prediction_list[sub])

        print("start convergence test")
        for loc_iter in range(30):
            for sub in range(n_sub):
                x_temp = []
                for i in range(history_level):
                    null_list = []
                    if(loc_iter < history_level - 1):
                        if(i < history_level - 1):    # first time 0,1, second time 0
                            for local_input in local_list[sub]:
                                null_list.append(dd_pod_coeffs_all[local_input][step + i])
                        else:                                # first time 2, second time 1,2, third time 0,1,2 and all the same after
                            for local_input in local_list[sub]:
                                null_list.append(prediction_update[local_input][-(history_level - i)])   # -3, -2, -1
                    else:
                        for local_input in local_list[sub]:
                            null_list.append(prediction_update[local_input][-(history_level - i)])
                    # for local_input in local_list[sub]:
                    #     # null_list.append(dd_pod_coeffs_all[local_input][step + i])
                    #     null_list.append(dd_pod_coeffs_all[local_input][step + i])
                    x_temp.append(null_list)
                x = Variable(torch.from_numpy(np.array(x_temp)).float())
    	        y_temp = []
                y_temp.append(dd_pod_coeffs_all[sub][step + history_level])
                y = Variable(torch.from_numpy(np.array(y_temp)).float())
                #print("x", x.size(), nodes_number*len(local_list[sub]))
                x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
                y = y.permute(-1,0).view(-1, nodes_number, 1)

                prediction = lstmNN[sub](x)
                prediction_update[sub].append(prediction.data.view(nodes_number).numpy())
                print("sub=", sub, np.max(prediction_update[sub][-1]-prediction_update[sub][-2]))
                #loss = loss_func(prediction, y)
                #optimizer[sub].zero_grad()
                #loss.backward()
                #optimizer[sub].step()
        prediction_list = []
        for sub in range(n_sub):
            x_np = []
            y_np = []
            x = []
            y = []
            for i in range(history_level):
                null_list = []
                for local_input in local_list[sub]:
                    null_list.append(prediction_update[local_input][-history_level+i])
                    #null_list.append(dd_pod_coeffs_all[local_input][step + i])
                # null_list.append(dd_pod_coeffs_all[local_list[sub][0]][step + i])
                x_np.append(null_list)
            x = Variable(torch.from_numpy(np.array(x_np)).float())
            y_np.append(dd_pod_coeffs_all[sub][step + history_level])
            y = Variable(torch.from_numpy(np.array(y_np)).float())
            x = x.permute(1,2,0).view(-1, nodes_number*len(local_list[sub]), history_level)
            # x = x.permute(1,2,0).view(-1, nodes_number, history_level)
            y = y.permute(-1,0).view(-1, nodes_number, 1)

            # x = x.permute(1,0).view(1, nodes_number*len(local_list[sub]), history_level)
            # y = y.permute(1,0).view(1, nodes_number, 1)
            prediction = lstmNN[sub](x)
            prediction_list.append(prediction.data.view(nodes_number).numpy())
            loss = loss_func(prediction, y)
            optimizer[sub].zero_grad()
            loss.backward()
            optimizer[sub].step()
    return lstmNN
