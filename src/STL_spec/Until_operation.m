function cell_until = Until_operation(cell_1, cell_2, start, End)
%%%% Until finction like every other function can be written based on 'and' and 'or' operations. This function shows how is this relation
%%%% phi1 U_[a,b] phi2    \exists t2\in[a,b] phi2(t2)==true & G_t\in[0,t2] phi1(t)==true
if start>0
    Cell = and_operation(  F_operation(   cell_2 , start,start  )  , G_operation(   cell_1  ,  0,start-1)   );
    for i = start+1:End
        Cell = or_operation(   Cell  , and_operation(   F_operation(   cell_2, i,i   )  ,  G_operation(   cell_1, 0,i-1 )  )  ) ;
    end
    cell_until = Cell;
else
    error('Lower time bound for Until operation can not be zero. ')
end
end
