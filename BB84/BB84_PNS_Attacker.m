function photon_nums = BB84_PNS_Attacker(origin_photon_nums)
    photon_nums = origin_photon_nums;
    L = length(origin_photon_nums);
    for i = 1:L
        if origin_photon_nums(i)<=1
            photon_nums(i) = 0;% keep the single photon or no photon
        else
            photon_nums(i) = origin_photon_nums(i)-1;% keep one photon
        end
    end
end