function H2 = Binary_Shannon_Entropy(x)
    if x==0||x==1
        H2 = 0;
    else
        H2 = -log2(x)*x-(1-x)*log2(1-x);
    end
end