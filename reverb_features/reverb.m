[r, fs] = generate_reverb('clean.wav', 'IR_new/h046_SuburbanGarage_4txts.wav');
audiowrite('ReverbAudio/SuburbanGarage_reverb.wav', r, fs);
soundsc(audioread('ReverbAudio/SuburbanGarage_reverb.wav'),fs);


