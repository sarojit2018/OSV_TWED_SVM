[y, Fs]=audioread("Vowel9.wav");
si = size(y);
speech=y(5000:40000);
window=50*0.001*Fs;
shift=round(window/4);
starter=1;
counter=0;
while counter<6
    short_signal = speech(starter:starter+window);
    starter = starter+shift;
    m=length(short_signal);
    n=pow2(nextpow2(m));
    spectrum = fft(short_signal, m);
    n=length(spectrum);
    f=(0:n-1) * (Fs/n);
    subplot(6, 2, 2*(6-counter)-1);
    plot(f(1:floor(n/9)), abs(spectrum(1:floor(n/9))));
    ceps=cceps(short_signal);
    n=length(ceps);
    t=(0:n-1)/(n*Fs);
    subplot(6, 2, 2*(6-counter));
    plot(t(1:floor(n/2)), ceps(1:floor(n/2)));
    counter=counter+1;
end
