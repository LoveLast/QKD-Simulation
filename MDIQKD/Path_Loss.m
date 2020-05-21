function photon_nums = Path_Loss(ori_nums,l,alpha0)
    t = 10^(-alpha0*l/10);% ¹âÏËÍ¸¹ıÂÊ
    photon_nums = ori_nums;
    for i = 1:length(ori_nums)
        photon_nums(i) = sum(randsrc(1,ori_nums(i),[0,1;1-t,t]));
    end
end