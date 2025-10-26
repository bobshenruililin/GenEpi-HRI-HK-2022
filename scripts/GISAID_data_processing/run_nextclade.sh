/opt/homebrew/bin/wget -O ~/software/nextclade/nextclade-x86_64-apple-darwin https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-apple-darwin

chmod 755 ~/software/nextclade/nextclade-x86_64-apple-darwin

~/software/nextclade/nextclade-x86_64-apple-darwin dataset get --name 'sars-cov-2' --output-dir $HOME'/software/nextclade/data/sars-cov-2'

~/software/nextclade/nextclade-x86_64-apple-darwin run \
	--in-order \
   	--input-dataset $HOME/software/nextclade/data/sars-cov-2 \
   	--output-all=$2 \
   	$1