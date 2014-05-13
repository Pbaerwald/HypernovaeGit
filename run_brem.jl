min_freq_ex = 15.3;
max_freq_ex = 22.0; 
num_freq = 500;
freq_ex = linspace(min_freq_ex,max_freq_ex,num_freq); 
freq = 10.^freq_ex; 

min_temp_ex = 7.0; 
max_temp_ex = 12.0; 
num_temp = 200; 
temp_ex = linspace(min_temp_ex,max_temp_ex,num_temp); 
temp = 10.^temp_ex; 

num_to_write = 10; 
results = cell(num_to_write); 

fp = open("parallel_bremsstrahlung_emissivity.dat","w");
index = 1; 
tic()
for i in 1:num_temp
	println("currently running temperature ",temp[i]); 
	index = i%num_to_write;
	if(index == 0)
		index = num_to_write; 
	end
	this_temp = temp[i]; 
	results[index] = pmap(c->specific_emissivity(c,this_temp),freq);
	
	if(index == num_to_write)
		println("writting..."); 
		for j in 1:num_to_write
			print(fp,results[j]'); 
		end
	end
end
 
for k in 1:index
	print(fp,results[k]'); 
end
toc(); 

close(fp); 
