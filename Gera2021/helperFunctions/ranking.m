function ranks=ranking(v) %determines the rank of the element in the vector;
[~,order]=sort(-v);
ranks=zeros(size(v));
ranks(order)=[1:max(size(v))];
end