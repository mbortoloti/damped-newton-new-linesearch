function [k] = snewtontest2(V,Gphip,P,metric,theta)
    indot = metric(Gphip,V,P);
    ng = sqrt(metric(Gphip,Gphip,P));
    nv = sqrt(metric(V,V,P));
    k = 1;
    if((indot + theta*ng*nv) > 0)
        k = 0;
    end
end
