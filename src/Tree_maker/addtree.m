function A = addtree(B, id, C)
n = nnodes(C);
ID = id;
for i = 1:n
    [B,new_id(i)] = addnode(B,ID,get(C,i));
    if i<n
        ID_C = getparent(C,i+1);
        ID = new_id(ID_C);
    end
end
    
A = B;

end
