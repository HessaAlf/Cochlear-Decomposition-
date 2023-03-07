function[angle, detail, approx] = CochlearDecomposition(signal, delta_theta) 
%% Cochlear Deccomposition. 


%Please cite H. Alfalahi, A. Khandoker, G. Alhussein, L. Hadjileontiadis " 
%Cochlear Decomposition: A Novel Bio-Inspired Multiscale Analysis 
%Framework" International Conference on Acoustics, Speech and Signal
%Processing, 2023. 


if nargin < 2
    
    delta_theta = 45; % default angular separation is 45 degree 
    
end 


full_angular_freq_range=20*(2*pi):2*pi:20000*(2*pi); %%% approximate hearing 
%frequency range in the human cochlea ~ 20 Hz - 20 kHz. 

theta = 0: delta_theta: 990;  


% position - frequency mapping 



for j=1:1:length(theta)
    
    angle(j) = theta(j);
    center_Freq(j)=165.4 .* ( ((3251.^(2.1)).* ((angle(j)+177.3).^(-2.1.*...
        1.149))) - 0.88);
    center_Freq(j) = 2*pi*center_Freq(j);
    tau(j)=(0.02)*center_Freq(j) ; 
    %[filter(j,:),comp(j,:),bm(j,:)]=GammaTone(x, center_Freq(j), fs);
end 



% create cochlear wavelets 


  for j=1:1:length(center_Freq)-1

    for i = 1:1:length(full_angular_freq_range)
   % for full_angular_freq_range= center_Freq(j):1:center_Freq(j+1)
   % or loop over angles, then determine the corresponding frequency
        %%these cutoffs sh[ould come from the angles
        %%what follows is that for every angle we have a scaling and
        %%wavelet function; which define the low pass and high pass
        %%quadrature mirror filters, respectively.
        if full_angular_freq_range(i) <= center_Freq(j) - tau(j)
            scaling(j,i)=1;
            
        elseif (center_Freq(j)-tau(j)< full_angular_freq_range(i)) && (full_angular_freq_range(i)<= center_Freq(j) + tau(j)) 
                scaling(j,i)= abs(cos(pi/2 * (beta((1/(2*tau(j)))*(full_angular_freq_range(i)-center_Freq(j)+tau(j))))));
                
                
        else
                scaling(j,i) = 0; 
                
        end
        
        if (center_Freq(j)-tau(j)<=full_angular_freq_range(i)) && (full_angular_freq_range(i)<= center_Freq(j) + tau(j) )
                wavelet(j,i)= abs(sin(pi/2 * (beta((1/(2*tau(j)))*(full_angular_freq_range(i)-center_Freq(j)+tau(j))))));
        
        
        elseif   center_Freq(j) + tau(j) <= full_angular_freq_range(i) && full_angular_freq_range(i)  <=  center_Freq(j+1) - tau(j+1)
            wavelet(j,i) = 1;
           
       elseif (center_Freq(j+1)-tau(j+1)<=full_angular_freq_range(i)) && (full_angular_freq_range(i) < center_Freq(j+1) + tau(j+1))
            
            wavelet(j,i)= abs(cos(pi/2 * (beta((1/(2*tau(j+1)))*(full_angular_freq_range(i)-center_Freq(j+1)+tau(j+1))))));
           
        
        else
            wavelet(j,i)=0;
        end
       
    end
  end
  
x=size(wavelet);   
  
for i=1:1:x(1)
    
wavelet_time(i,:) =irfft(wavelet(i,:),length(wavelet(i,:)));
N=length(wavelet_time(i,:)); 
wavelet_time=wavelet_time(i,:);
%try checking if N is odd or even.
wavelet_time1(i,:)=wavelet_time(1:(N-1)/2);
wavelet_time2(i,:)=wavelet_time((N)/2:N);

wavelet_full(i,:) = [wavelet_time2(i,:), wavelet_time1(i,:)]; 
%wavelet_full(i,:)=wavelet_full(i,:)/max(wavelet_full(i,:));
end 


  
for i=1:1:x(1)
    
scaling_time(i,:) =irfft(scaling(i,:),length(scaling (i,:)));
 N=length(scaling_time(i,:)); 
 scaling_time=scaling_time(i,:);
% %try checking if N is odd or even.
 scaling_time1(i,:)=scaling_time(1:(N-1)/2);
 scaling_time2(i,:)=scaling_time((N)/2:N);
% 
 scaling_full(i,:) = [scaling_time2(i,:), scaling_time1(i,:)]; 
 %scaling_full(i,:)=scaling_full(i,:)/max(scaling_full(i,:));
end 



input=signal;

for i=x(1):-1:1
    
approx(i,:)=conv(input,scaling_full(i,:), 'same');

detail(i,:)=conv(input,wavelet_full(i,:), 'same'); 

%input=approx(i,:); 
end 



end 