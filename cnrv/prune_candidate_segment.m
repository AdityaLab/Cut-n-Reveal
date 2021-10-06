function [f] = prune_candidate_segment(segLenThresh,segmentIndices,numCluster,ttlTime,isV)

f = 1;
sz = size(segmentIndices);
if (sz+1)<numCluster
    f= 0;
else
    if (isV)
        prevT= segmentIndices(1);
        for t=2:size(segmentIndices,2)
            diff = segmentIndices(t)-prevT;
            prevT = segmentIndices(t);
            if diff<(ttlTime*segLenThresh)
                f =0;
                break;
            end
        end
    end
end
end
