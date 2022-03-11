import deepmd.DeepDipolestone     as DP
import deepmd.DeepPrimitiveMoment as QP
import numpy as np

dp = DP('graph_m2.pb')
qp = QP('graph_m4.pb') 

num_water = 64
num_atoms = num_water*3

L = 13.023575
cell = np.diag(L*np.ones(3)).reshape([1,-1])
atype = [0, 1, 1]*int(num_water)

frame0 = 1
frame1 = 30000

file_dip = open("m2.1","a")
file_qua = open("m4.1","a")

with open('water.xyz') as file_xyz:

    ndx_skip = 1
    for ndx_frame in range(1,frame0):
        head_line  = (ndx_frame-1)*(num_atoms+2)+1
        tail_line  = head_line + num_atoms+2 - 1
        for ndx_line in range(head_line,tail_line+1):
            line = file_xyz.readline()
            ndx_skip = ndx_skip+1

    print('ndx_skip = ', ndx_skip)


    for ndx_frame in range(frame0,frame1+1):

        coord_list = []
        head_line  = (ndx_frame-frame0)*(num_atoms+2)+ndx_skip
        tail_line  = head_line + num_atoms+2 - 1

        ndx_record = 1
        for ndx_line in range(head_line,tail_line+1):
            line = file_xyz.readline()
            if 2 < ndx_record and ndx_record <= num_atoms+2:
                xyz = []
                cont = line.split()
                for u in cont[1:4]:
                    trimed = u.strip()
                    xyz.append(float(trimed))
                coord_list.append(xyz)
            ndx_record = ndx_record + 1

        coord = np.array(coord_list).reshape([1,-1])
        dip = dp.eval(coord,cell,atype)
        qua = qp.eval(coord,cell,atype)
        for i in range(num_water):
            print(ndx_frame,i+1, dip[0][i][0],dip[0][i][1],dip[0][i][2])
            file_dip.write('%d %d %8.4f %8.4f %8.4f\n'%(ndx_frame, i+1, dip[0][i][0],dip[0][i][1],dip[0][i][2]))
            print(ndx_frame,i+1,qua[0][i][0],qua[0][i][1],qua[0][i][2])
            print(ndx_frame,i+1,qua[0][i][3],qua[0][i][4],qua[0][i][5])
            print(ndx_frame,i+1,qua[0][i][6],qua[0][i][7],qua[0][i][8])
            file_qua.write('%d %d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' \
                    %(ndx_frame, i+1, qua[0][i][0],qua[0][i][1],qua[0][i][2], \
                                      qua[0][i][3],qua[0][i][4],qua[0][i][5], \
                                      qua[0][i][6],qua[0][i][7],qua[0][i][8]))

file_dip.close()
file_qua.close()



