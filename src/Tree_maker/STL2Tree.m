function Tree = STL2Tree(str)

%%%% please avoid using unnecessary paranthesis it results in deep trees %%%%

%%% In this function we separate the STL formula based on the paranthesis. for example:

%%% ((a|c) & d & (e|f))  ----->   1- (a|c)   2- d  3- (e|f)   ---> It also put these 3 terms on a logic tree.             &       <------   root 
%%%                                                                                                                 /     \   \
%%%                                                                                                               (a|c)    d  (e|f)   <----  childeren
%%% in one attempt it splits the terms in paranthesis. In the next attempt it goes again through the splitted parts and check if they are again splittable 
%%% with the paranthesis or not. If YES then it serarates the string again
%%% and creates a new Tree for that, and replace that tree with the splitted string. On the main Tree. We continue until no paranthesis exists any more.
%%%  Therefore         &    <--- root
%%%               /     \   \
%%%              |      d    |   <----  childeren
%%%             / \         / \
%%%            a   c       e   f  <---- childeren


Tree =  separate_operator(str);

dd = true;
while dd  %%% In this while loop we perform the steps introduced in comments --->  Separation until no separation is possible. 
    n = nnodes(Tree);
    for i = 1:n
        str = get(Tree,i);
        if length(str)>1
            if str(1) =='('
                Tree1 = separate_operator(str(2:end-1));
                Tree =  replace_node_with_tree(Tree, i, Tree1);
            else
                Tree1 = separate_operator(str);
                Tree =  replace_node_with_tree(Tree, i, Tree1);
            end
        
        end
    end
    m = nnodes(Tree);
    if m==n
        dd = false;
    end
end

%%% remove unnecessary 'and' and 'or' operations
%%%  After creation of Tree we attempt to apply logical filters, to decrease the depth of the tree if it is possible.

n1 = nnodes(Tree);
n = n1;
i = 0;
while i<n    %%% The first filter, checks there should not be a parant of form '&' that has childeren of form '&'. If it finds it, then it removes the lower and and bring up its childeren. The same thing happens for '|'.
    i = i+1;
    if strcmp(get(Tree,i),'and')
        kids = getchildren(Tree,i);
        m1 = length(kids);
        m = m1;
        j = 0;
        while j<m
            j = j+1;
            if strcmp(get(Tree,kids(j)),'and')
                Parent = getparent(Tree,kids(j));
                Tree = removenode(Tree,kids(j));
                kids = getchildren(Tree, Parent);
                m = length(kids);
                n = nnodes(Tree);
                j = 0;
                i = 0;
            end
        end
    end
end
n1 = nnodes(Tree);
n = n1;
i = 0;
while i<n
    i = i+1;
    if strcmp(get(Tree,i),'or')
        kids = getchildren(Tree,i);
        m1 = length(kids);
        m = m1;
        j = 0;
        while j<m
            j = j+1;
            if strcmp(get(Tree,kids(j)),'or')
                Parent = getparent(Tree,kids(j));
                Tree = removenode(Tree,kids(j));
                kids = getchildren(Tree, Parent);
                m = length(kids);
                n = nnodes(Tree);
                j = 0;
                i = 0;
            end
        end
    end
end

%%% The second filter check the siblings and removes repeatative sibling leaves.
leave = true;
while leave
    IDs  =  findleaves(Tree);
    n = length(IDs);
    leave = false;
    for i = 1:n-1
        for j = i+1:n
            if strcmp(get(Tree,IDs(i)),get(Tree,IDs(j)))  && ismember(  IDs(j)   ,   getsiblings(Tree,IDs(i))    )
                leave =  true;
            end
            if leave
                break;
            end
        end
        if leave
            break;
        end
    end
    
    if leave
        Tree = removenode(Tree,IDs(j));
    end
end
            


decision = nnodes(Tree);
if decision <500
    figure
    plot(Tree)
end

% n=nnodes(Tree);
% for i=1:n
%     if isleaf(Tree,i)
%         Tree=set(Tree,i,str2num(get(Tree,i)));
%     end
% end
    
end
    