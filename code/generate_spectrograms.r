window <- 256
run_length <- 4096
create_spectrogram <- function(brain, brain_filename, electrode) {
    cat(paste('Processing data for spectrogram:',brain_filename,'\n'))

    # Load the electrode series.
    brain_series <- brain[1:(run_length+window),electrode+1]
    
    # Center the brainwaves using the moving average.
    cat(paste('Centering brainwaves\n'))
    brain_series_avg <- c()
    for(i in 1:length(brain_series)) {
        brain_series_avg <- c(brain_series_avg, mean(brain_series[i:i+10],na.rm=TRUE))
    }
    brain_series_center <- brain_series-brain_series_avg

    # Fourier transform of the brainwaves.
    library(stats)
    freqs <- c()
    cat(paste('Getting Fourier transform\n'))
    for(i in 1:run_length) {
        # Get log of frequency amplitudes.
        amps <- log(Mod(fft(brain_series_center[i:(i+window-1)])))
        # Center frequency amplitudes.
        amps_center <- (amps-min(amps))/(max(amps)-min(amps))
        freqs <- rbind(freqs, amps_center)
    }

    brain_path <- paste('../images/',brain_filename,sep='')
    png(brain_path)
    b <- window/2+1 # Half of the spectrum is redundant
    cat(paste('Drawing spectrogram in', brain_path, '\n'))
    image(freqs[,1:b],col=rainbow(run_length),xaxt='n',yaxt='n')
    tmp <- dev.off()
}

# Load brainwave datasets
filename <- "../data/brain_tapping.csv" # Tapping
brain_tapping <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))
filename <- "../data/brain_inactive.csv" # Inactive
brain_inactive <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))
filename <- "../data/brain_tetris_slow.csv" # Slow Tetris
brain_tetris_slow <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))
filename <- "../data/brain_tetris_fast.csv" # Fast Tetris 
brain_tetris_fast <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))
filename <- "../data/brain_furrow.csv" # Furrowed brow
brain_furrow <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))
filename <- "../data/brain_tension.csv" # Face and jaw tensed
brain_tension <- read.csv(filename, comment.char="%", stringsAsFactors=FALSE)
cat(paste('loaded',filename,'\n'))

# Create the spectrograms.
create_spectrogram(brain_tapping,     "brain_tapping.png",     1)
create_spectrogram(brain_inactive,    "brain_inactive.png",    1)
create_spectrogram(brain_tetris_slow, "brain_tetris_slow.png", 1)
create_spectrogram(brain_tetris_fast, "brain_tetris_fast.png", 1)
create_spectrogram(brain_furrow,      "brain_furrow.png",      1)
create_spectrogram(brain_tension,     "brain_tension.png",     1)