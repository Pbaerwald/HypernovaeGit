fp = open("test.txt","r");
result = readline(fp);
println(result);
println(typeof(result));

#while !eof(fp)
#	append!(readline(fp),result); 
#end
close(fp);

for i in 1:length(result)
	print(result[i]); 
end
