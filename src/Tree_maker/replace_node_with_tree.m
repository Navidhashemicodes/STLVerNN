function A = replace_node_with_tree(B, id, C)

%%%% This function Takes the tree (B) and selects the node with position(id). Then it takes the new tree C and replace it with node in position (id). Finally we return the new tree as A.

n = nnodes(C);
B = set(B,id,get(C,1));
ID = id;
new_id(1) = id;
for i = 2:n
    [B,new_id(i)] = addnode(B,ID,get(C,i));
    if i<n
        ID_C = getparent(C,i+1);
        ID = new_id(ID_C);
    end
end
    
A = B;

end