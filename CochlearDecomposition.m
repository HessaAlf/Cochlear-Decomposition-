function[angle, detail, approx] = CochlearDecomposition(signal, delta_theta) 
%% Cochlear Deccomposition. 


%Please cite H. Alfalahi, A. Khandoker, G. Alhussein, L. Hadjileontiadis " 
%Cochlear Decomposition: A Novel Bio-Inspired Multiscale Analysis 
%Framework" International Conference on Acoustics, Speech and Signal
%Processing, 2023. 


if nargin < 2
    
    delta_theta = 45; % default angular separation is 45 degree 
    
end 


full_angular_freq_range=10*(2*pi):2*pi:20000*(2*pi); %%% approximate hearing 
%frequency range in the human cochlea ~ 20 Hz - 20 kHz. 

theta = 0: delta_theta: 990;  


% position - frequency mapping 



for j=1:1:length(theta)
    
    angle(j) = theta(j);
    center_Freq(j)=165.4 .* ( ((3251.^(2.1)).* ((angle(j)+177.3).^(-2.1.*...
        1.149))) - 0.88);
    center_Freq(j) = 2*pi*center_Freq(j);
    if j ~= 1
    x(j) = (center_Freq(j)-center_Freq(j-1))/(center_Freq(j)+center_Freq(j-1)); 
    end
    tau(j)=(0.1)*center_Freq(j); 
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
                scaling(j,i)= abs(cos(pi/2 * (beta((1/(2*tau(j)))*(full_angular_freq_range(i)-(center_Freq(j))+tau(j))))));
                
                
        else
                scaling(j,i) = 0; 
                
        end
        
        if (center_Freq(j)-tau(j)<=full_angular_freq_range(i)) && (full_angular_freq_range(i)<= center_Freq(j) + tau(j))
                wavelet(j,i)= abs(sin(pi/2 * (beta((1/(2*tau(j)))*(full_angular_freq_range(i)-(center_Freq(j))+tau(j))))));
        
        
        elseif   center_Freq(j) + tau(j) <= full_angular_freq_range(i) && full_angular_freq_range(i)  <=  center_Freq(j+1) - tau(j+1)
            wavelet(j,i) = 1;
           
       elseif (center_Freq(j+1)-tau(j+1)<=full_angular_freq_range(i)) && (full_angular_freq_range(i) < center_Freq(j+1) + tau(j+1))
            
            wavelet(j,i)= abs(cos(pi/2 * (beta((1/(2*tau(j+1)))*(full_angular_freq_range(i)-(center_Freq(j+1))+tau(j+1))))));
           
        
        else
            wavelet(j,i)=0;
        end
       
    end
  end
  
x=size(wavelet);   
  
for i=1:1:x(1)
 Y1 = wavelet(i,:);
 %Y2= [Y1(1) Y1(2:end)/2 fliplr(conj(Y1(2:end)))/2];
 Y2= [Y1(1) Y1(2:end)/2];
 %normalization_factor = max(abs(wavelet(i,:)));
% magnitude_spectrum = abs(normalization_factor); 

wavelet_time(i,:) =real(ifftshift(ifft(Y2,'symmetric')));
N=length(wavelet_time(i,:)); 
%wavelet_time(i,:) =conj(wavelet_time(i,:)); 
%wavelet_times=wavelet_time(i,:);
%try checking if N is odd or even.
% wavelet_time1=wavelet_times(1:(N-1)/2);
% wavelet_time2=wavelet_times((N)/2:end);
% 
% wavelet_full(i,:) = [wavelet_time2, wavelet_time1]; 
% wavelet_full(i,:)=wavelet_full(i,:)/max(wavelet_full(i,:));
end 


  
for i=1:1:x(1)
 Y1 = scaling(i,:);
 Y2= [Y1(1) Y1(2:end)/2 fliplr(conj(Y1(2:end)))/2];   
scaling_time(i,:) =real(ifftshift(ifft(Y2)));
%  N=length(scaling_time(i,:)); 
%  scaling_time(i,:) = N*scaling_time(i,:); 
%  scaling_times=scaling_time(i,:);
% % %try checking if N is odd or even.
%  scaling_time1(i,:)=scaling_times(1:(N-1)/2);
%  scaling_time2(i,:)=scaling_times((N)/2:N);
% % 
%  scaling_full(i,:) = [scaling_time2(i,:), scaling_time1(i,:)]; 
%  %scaling_full(i,:)=scaling_full(i,:)/max(scaling_full(i,:));
end 

% perhaps we try with the complex wavelet 


input=signal;

%energy = norm(wavelet_time(1,:)/max(wavelet_time(1,:)));

for i=x(1):-1:1
    
    approx(i,:)=conv(input,scaling_time(i,:), 'same');

    detail(i,:)= conv(input,wavelet_time(i,:), 'same'); 

%input=approx(i,:); 
end 




end 
