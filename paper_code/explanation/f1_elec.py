import sys

def f1_elec(seg_dir, segfile, gt_file, error, save_file):
    #get the segmentation, the ground truth cut, the tolerance, calculate the f1 score. 
    print 'reading ground truth'
    with open(seg_dir + gt_file) as f:
        line = f.readline()
    gt_c = [int(x) for x in line.strip().split(',')]
    print 'reading segmentations'
    with open(seg_dir + segfile) as f:
        line = f.readline()
    S = [int(x) for x in line.strip().split(',')]
    print 'calculate f1'
    tp = 0.0
    fp = 0.0
    fn = 0.0
    print("Number of Predicted Segments",len(S),"Ground Truth Segments",len(gt_c))
    for s in S:
        if hit(s, gt_c, error):
            tp += 1
        else:
            fp += 1
    for c in gt_c:
        if not hit(c, S, error):
            fn += 1
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = 2 * (prec * rec) / (prec + rec)
    print 'save results'
    print "f1 score",f1,"prec",prec,"rec",rec,"tp",tp,"fp",fp
    sf = open(seg_dir + save_file, 'wb')
    sf.write(str(f1))
    sf.close()

def hit(s, li, e):
    #check if s is in the set l with tolerance e
    for l in li:
        if abs(l - s) <= e:
            return True
    return False

if __name__ == '__main__':
    #Pass in lambda_1,lambda_2,numcluster,seg_dir
    seg_dir = str(sys.argv[5])
    l1 = str(sys.argv[1])
    l2 = str(sys.argv[2])
    l3 = str(sys.argv[3])
    numcluster=str(sys.argv[4])
    error = int(sys.argv[7]) 
    gt_file = str(sys.argv[6])
    #gt_file = 'gt_cuts.txt'#ground truth cut point file
    segfile = 'segV_lam1_' + l1 + '_lam2_' + l2 + '_lam3_' + l3 +'_clusV_'+numcluster+'.csv' 
    #segfile = str(sys.argv[8])
    save_file = 'f1_' + l1 + '_' + l2 + '_' + l3+ '.txt'
    f1_elec(seg_dir, segfile, gt_file, error, save_file)
