import sys
import glob
import os

def f1_elec(segfile, gt_file, error, save_file):
    #get the segmentation, the ground truth cut, the tolerance, calculate the f1 score. 
    #print 'reading ground truth'
    with open(gt_file) as f:
        line = f.readline()
    gt_c = [int(x) for x in line.strip().split(',')]
    #print 'reading segmentations'
    with open(segfile) as f:
        line = f.readline()
    S = [int(x) for x in line.strip().split(',')]
    print 'calculate f1'
    tp = 0.0
    fp = 0.0
    fn = 0.0
    seen=list()
    seenfps=list()
    for idx,s in enumerate(S):
        _idx,ishit=hit(s,gt_c,error)
        if ishit and (_idx not in seen):
            tp += 1
            seen.append(_idx)
        else:
            if _idx not in seenfps:
                seenfps.append(_idx)
                fp += 1

    fp=abs(len(S) - tp) 
    for c in gt_c:
        _idx,ishit=hit(c,S,error)
        if not ishit:
            fn += 1
    fn=max(0,len(gt_c)-tp)
    prec = tp / max(1,(tp + fp))
    rec = tp / max(1,(tp + fn))
    f1 = (2 * (prec * rec)) / max(1,(prec + rec))
    #print 'save results'
    #if f1>-1:
    #print("Number of Predicted Segments"+str(len(S))+"Ground Truth Segments"+tr(len(gt_c)))
    print "f1 score:",f1," prec:",prec," rec:",rec," tp:",tp," fp:",fp," fn",fn
    sf = open(seg_dir + save_file, 'wb')
    sf.write(str(f1))
    sf.close()

def hit(s, li, e):
    #check if s is in the set l with tolerance e
    for l in li:
        if abs(l - s) <= e:
            return l,True
    return None,False

if __name__ == '__main__':
    #Pass in lambda_1,lambda_2,numcluster,seg_dir
    seg_dir = '../result/cnrV/wjr1/'#str(sys.argv[1])
    gt_file = '../result/cnrV/wjr1_gt.txt'#ground truth cut point file
    for _file in glob.glob(seg_dir+'osc_segment_indices_*.csv'):
        #l1 = str(sys.argv[1])
        #try:
            item=os.path.basename(_file)
            l1=item.split('_')[5]
            l2=item.split('_')[8]
            #print("Tokens",item.split('_'))
            numcluster=item.split('_')[10]
            numiter=300#item.split('_')[10]
            #l2 = str(sys.argv[2])
            #numcluster=str(sys.argv[3])
            #numiter=str(sys.argv[5])
            #error=79 #chicken19
            #error=36 #nilm1
            #error=16 #chicken21
            error=50 #wjr1,grandmal,synthetic
            #error = 15  #wjr2
			#error = 10 #Window for Chicken
            #error = 5 #Window for Walk Jog Run
            #error = 22 #Window for NILM
            #gt_file = 'gt_cuts.txt'#ground truth cut point file
            #segfile = 'osc_segment_indices_lambda_1_' + l1 + '_lambda_2_' + l2 + '_numiter_'+numiter+'_numcluster_'+numcluster+'.csv' 
            save_file = 'f1_' + l1 + '_' + l2 + '.txt'
            print "L1",l1,"L2",l2,"NumIter",numiter,"Num Cluster",numcluster
            f1_elec(_file, gt_file, error, save_file)
            print("\n=============\n")
        #except:
         #   continue
