#!/usr/bin/perl
#INPUT: ldel/ltot treatment replicate#
#input example: perl multideme.pl 0.8 1 0 0: high 80% initial Ldel with changing environment and a unimodal environment fitness function
use strict;
use Math::Random qw(random_gamma random_normal random_multinomial);

my $num_del_loci =500;
my $gamma = 20;#0.007* $num_del_loci / 0.5;
my $num_trait_loci = 10;
my $num_traits = 1;	
#chomp $row;
						

my $num_demes = 1;
my $L_equilibrium = 0.8;
my $pop_magnitude = 4;#log(500)/log(10);
my $Ne = 10**$pop_magnitude;
my $mu = 10**-8;
my $mu_a= (1 + 7/9 + 7/9) * $mu;
my $delta = 10**-2.5;
my $sigma_f = 0.5;
my $sigmam = 2;
my $sigma_e = 4;
my $p_del = 0.4;
my $p_minus_del = 0.1;
my $environment_duration = 15000;
my $generations = 1000000;# $environment_duration*10000;
my @Ldel_init = ($ARGV[0] * $num_del_loci)x $num_demes;
my $replicate = $ARGV[3];
my $environmental_change = 0;#$ARGV[1];
my $burnin = 0;#$environment_duration*5000;
my $num_printouts =1000;
my $mode = $ARGV[1];
my $env_bias = 1/5;
my $trait_mut_bias = 1/50;
my $using_mode;
my $cooption = 0;#$ARGV[2];
my $migration_prob = 0.0;
my $record_rho = 1;
my @high_rho_values;
my $testing_capacitance = 0;
my $bumpiness = $ARGV[2];

if ($ARGV[1] == 0){
	$using_mode = "none";
}elsif ($ARGV[1] == 1){
	$using_mode = "env_ch";$environmental_change = 1;
}elsif ($ARGV[1] == 2){
	$using_mode = "both";$testing_capacitance = 1;
}elsif ($ARGV[1] ==3){	
	$using_mode = "low";$testing_capacitance = 1;
}elsif ($ARGV[1]==4){
	$using_mode = "high";$testing_capacitance = 1;
}elsif ($ARGV[1]==5){
	$using_mode = "time";$testing_capacitance = 1;
}elsif ($ARGV[1]==6){
	$using_mode = "staticLdel_ch";$environmental_change = 1;
}elsif ($ARGV[1]==7){
	$using_mode = "staticLdel_noch";
}

$generations += $burnin;	
my @high_rho_values;

my @total_bouts;
my @total_high_rho;
if ($using_mode eq "time" || $using_mode eq "high" || $using_mode eq "both"){

	for(my $ldel = 400; $ldel > 0; $ldel -=25){
		my $ldel_index = int($ldel / 25);
		$total_bouts[$ldel_index] = 0;
		$total_high_rho[$ldel_index] = 0;
		for(my  $rep = 0; $rep <= 2; $rep++){
			my $high_rho_input = "rhooutput_4_" . $ldel . "_staticLdel_ch_" . $rep . ".txt";
			#print "opening file $high_rho_input\n";

			open (my $highrhofile, '<', $high_rho_input) or die "Could not open file '$high_rho_input' $!";
			my $current_bout = -1;
			while (my $row = <$highrhofile>) {
				chomp $row;
				my @values=split("\t",$row);
				if($values[0]!=$current_bout and int($values[1]/15000) > 6000){
					$current_bout++;
					$total_bouts[$ldel_index]++;
				}
				if(int($values[1]/15000) > 6000){
					$total_high_rho[$ldel_index]++;
				}
			}
		}
		#print $ldel . "\t" . $ldel_index  . "\t" . $total_bouts[$ldel_index] . "\n";
	}
}
foreach (my $i = 0; $i < 17 ;$i++){
	#print $i . "\t" . $total_bouts[$i] . "\t" . $total_high_rho[$i]. "\n";
}
if ($using_mode eq "time"){
	for(my  $ldel = 400; $ldel > 0; $ldel -=25){
		my $ldel_index = int($ldel / 25);
		my $i = 0;
		for(my  $rep = 0; $rep <= 2; $rep++){
			my $total_bout = -1;
			my $high_rho_input = "rhooutput_4_" . $ldel . "_staticLdel_ch_" . $rep . ".txt";
			

			open (my $highrhofile, '<', $high_rho_input) or die "Could not open file '$high_rho_input' $!";
			my $current_bout = -1;
			my $record_bout = 0;
			while (my $row = <$highrhofile>) {
				chomp $row;
				my @values=split("\t",$row);
				if($values[0]==$current_bout){
					$i++;
				}else{
					if(rand()<200/$total_bouts[$ldel_index] and int($values[1]/15000) > 6000){
						$record_bout=1;
						$total_bout++;
					}else{$record_bout=0;}
					$current_bout++;
					$i=0;
				}
				if ($record_bout==1){
					$high_rho_values[$ldel_index][$total_bout][$i]=$values[3];# print $ldel_index. "\t" . $total_bout . "\t" . $i . "\t" .$high_rho_values[$ldel_index][$total_bout][$i] . "\n";
				}
			}
		}
	}
}

if ($using_mode eq "high" || $using_mode eq "both"){
	for(my  $ldel = 400; $ldel > 0; $ldel -=25){
		my $ldel_index = int($ldel / 25);
		my $i = 0;
		for(my  $rep = 0; $rep <= 2; $rep++){
			my $high_rho_input = "rhooutput_4_" . $ldel . "_staticLdel_ch_" . $rep . ".txt";
			#print "opening file $high_rho_input\n";

			open (my $highrhofile, '<', $high_rho_input) or die "Could not open file '$high_rho_input' $!";	
			while (my $row = <$highrhofile>) {
				chomp $row;
				my @values=split("\t",$row);
				if(rand()<2000/$total_high_rho[$ldel_index] and int($values[1]/15000) > 6000){
					$high_rho_values[$ldel_index][$i]=$values[3];
					$i++;
				}
			}
		}
		#print $i . "\n";
	}
}
#print "done reading\n";
#sleep 5;
my $beta_mut_bias = 50;

for (my $holder = 0; $holder < 1; $holder ++){#(my $replicate = 0; $replicate < 2; $replicate++){#(my $pop_magnitude  = 4; $pop_magnitude <= 4; $pop_magnitude +=1){#
	my $filename  = "evol_" . $pop_magnitude . "_" . $Ldel_init[0] . "_" . $using_mode . "_" . $bumpiness . "_" . $replicate . ".txt";	


	open OUTPUT, '>', $filename or die "Can't open $filename - $!\n";
	printf OUTPUT "Ne\tLdel_init\tenv_ch\tgen\tLdel\tlog_rho\tenv_opt\talpha_sum\ttraitdel\tbeta_sum\tnum_highrho\n"; 
	if ($record_rho == 1){
		my $rhofilename  = "rhooutput_" . $pop_magnitude . "_" . $Ldel_init[0] . "_" . $using_mode . "_" . $replicate . ".txt";	
		open RHOOUTPUT, '>', $rhofilename or die "Can't open $rhofilename - $!\n";
	}

	my $rhofreqfilename  = "rhofreq_" . $pop_magnitude . "_" . $Ldel_init[0] . "_" . $using_mode . "_" . $replicate . ".txt";	
	open RHOFREQOUTPUT, '>', $rhofreqfilename or die "Can't open $rhofreqfilename - $!\n";

	$Ne = 10**$pop_magnitude;
	my $meanlogpfix = 0;
	my $sumldelmutants = 0;
	my $sumldelminusmutants = 0;
	my $sumldelminusfixations = 0;


	my $ldel_mut_freq = $num_del_loci * 30*10**-8;
	my $rho_mut_freq = 10**-5;
	my $trait_mut_freq = $num_trait_loci * 330 * 10**-8;
	
	if ($using_mode eq "staticLdel_ch"||$using_mode eq "staticLdel_noch"){
		$ldel_mut_freq = 0;
	}
	if ($using_mode eq "both" ||$using_mode eq "low"||$using_mode eq "time"||$using_mode eq "high"){
		$num_trait_loci = 1;
	}

	my $tot_mut_frequency = $ldel_mut_freq + $rho_mut_freq + $trait_mut_freq;
	my $ldel_mut_prob = $ldel_mut_freq / $tot_mut_frequency;#0.681818182
	my $rho_mut_prob = $rho_mut_freq/ $tot_mut_frequency;#0.045454545
	my $trait_mut_prob = $trait_mut_freq/ $tot_mut_frequency;#0.272727273	
	my @mut_probs = ($ldel_mut_prob, $rho_mut_prob, $trait_mut_prob);

	my @alpha = [];
	my @beta = [];
	my @genotype = [];
	my @traitdel = [];
	my @environmental_optimum = [];
	my @Ldel = [];
	my @rho = [];
	for (my $deme =0; $deme < $num_demes; $deme++){
		for (my $i =0; $i < $num_trait_loci; $i++){
			for (my $j =0; $j < $num_traits; $j++){
				$alpha[$deme][$i][$j]=0.0;#random_normal(1,0,2);
				$beta[$deme][$i][$j]=0.0;#random_normal(1,0,2);
			}
		}
		for (my $i = 0; $i < $num_trait_loci; $i++){
			if ($i < $ARGV[0] * $num_trait_loci){
				$traitdel[$deme][$i] = 1;
			}else{
				$traitdel[$deme][$i] = 0;
			}
		}
		for (my $i = 0; $i < $num_traits; $i++){
			$environmental_optimum[$deme][$i]=0;
		}
	}
		
	for (my $deme =0; $deme < $num_demes; $deme++){
		my $sum_trait_del = 0;
		for (my $i = 0; $i < $num_trait_loci; $i++){
			$sum_trait_del+=$traitdel[$deme][$i];
		}
		$Ldel[$deme]=$Ldel_init[$deme]-$sum_trait_del;
		if ($Ldel[$deme] < 0){$Ldel[$deme] = 0;}
	}

	for (my $deme =0; $deme < $num_demes; $deme++){
		my $max_fit = 0;
		my $best_rho = 0;
		my @curr_alpha;
		my @curr_beta;
		my @curr_traitdel;
		for (my $i =0; $i < $num_trait_loci; $i++){
			for (my $j =0; $j < $num_traits; $j++){
				$curr_alpha[$i][$j]=$alpha[$deme][$i][$j];
				$curr_beta[$i][$j]=$beta[$deme][$i][$j];
			}
			$curr_traitdel[$i] = $traitdel[$deme][$i];
		}
		my @curr_env_opt;
		for (my $i = 0; $i < $num_traits; $i++){
			$curr_env_opt[$i] = $environmental_optimum[$deme][$i];
		}
		for (my $i = 0.2; $i > 10**-13; $i *= 0.9){
			my @tempgenotype = ($Ldel[$deme], $i);
			my @curr_env_opt = 0;
		#	print $i . "\t";
			my $fitness = calculate_fitness(\@tempgenotype, \@curr_alpha, \@curr_beta, \@curr_traitdel, \@curr_env_opt, 0);
			if ($fitness > $max_fit){
				$best_rho = $i;
				$max_fit = $fitness;
				#print $fitness . "\n";

			}
		}
		$rho[$deme]=$best_rho;
	}
	my $rho_record = 0;
	my $rho_record_num = 0;
	my @rho_recorder;
	my @beta_recorder;
	my @gen_recorder;
	my @alpha_recorder;
	my @env_recorder;
	my @ldel_recorder;
	my @high_rho_duration;
	my $high_rho_bout = 0;	
	my @rho_frequency = (0) x 100;
	my $high_rho_sample_number = 0;
	my $high_rho_bout_number = 0;
	my $num_high_rho_values_used = 0;
	my $Ldel_index;
	for (my $gen = 0; $gen < $generations+1; $gen++){ 
		#print $gen ."\t";
		my @new_Ldel;
		my @new_rho;
		my @new_alpha;
		my @new_beta;
		my @new_traitdel;
		my @migration;
		my @env_opt;
		for (my $deme =0; $deme < $num_demes; $deme++){	
			$migration[$deme]=0;
			my $source_deme = $deme;
			if(rand()<$migration_prob){
				$source_deme=int(rand()*$num_demes);
				$migration[$deme]=1;
			}
			$new_Ldel[$deme] = $Ldel[$source_deme];
			$new_rho[$deme] = $rho[$source_deme];	
			for (my $i = 0; $i < $num_trait_loci; $i++){				
				for (my $j = 0; $j < $num_traits; $j++){
					$new_alpha[$deme][$i][$j] = $alpha[$source_deme][$i][$j];
					$new_beta[$deme][$i][$j] = $beta[$source_deme][$i][$j];
				}
				$new_traitdel[$deme][$i] = $traitdel[$source_deme][$i];
			}
		}

		for (my $deme =0; $deme < $num_demes; $deme++){
			my @mut_genotype = ($new_Ldel[$deme], $new_rho[$deme] );
			my @curr_alpha;
			my @curr_beta;
			my @curr_traitdel;
			my @mut_alpha = @curr_alpha;
			my @mut_beta = @curr_beta;
			my @mut_traitdel = @curr_traitdel;			

			for (my $i = 0; $i < $num_traits; $i++){
				$env_opt[$i]=$environmental_optimum[$deme][$i];
			}
			
			for (my $i = 0; $i < $num_trait_loci; $i++){
				for (my $j = 0; $j < $num_traits; $j++){
					$curr_alpha[$i][$j] = $new_alpha[$deme][$i][$j];
					$curr_beta[$i][$j] = $new_beta[$deme][$i][$j];
					$mut_alpha[$i][$j] = $new_alpha[$deme][$i][$j];
					$mut_beta[$i][$j] = $new_beta[$deme][$i][$j];
				}
				$curr_traitdel[$i] = $new_traitdel[$deme][$i];
				$mut_traitdel[$i] = $new_traitdel[$deme][$i];

			}			

			my $ldel_mutant = 0;
			my $trait_mut = 0;
			my $sum_trait_del;
			my $mutation;
			foreach(@mut_traitdel){$sum_trait_del+=$_;}
			if($migration[$deme] ==0){
				my @mut_rand = random_multinomial(1, @mut_probs);
				if ($mut_rand[0]==1){
					$mutation = "Ldel";
					if (rand() < $mut_genotype[0] / ($num_del_loci - $num_trait_loci)){
						if(rand()<$p_minus_del){
							$sumldelmutants++;
							$mut_genotype[0]--;
 							$ldel_mutant = 1; $mutation = "Ldel-";
 							$sumldelminusmutants++;
						}	
					}elsif (rand()<$p_del){
						$mut_genotype[0]++;
						$ldel_mutant = 1; $mutation = "Ldel+";
						$sumldelmutants++;
					}
				}
				if ($mut_rand[1]==1){
					$mutation = "rho";
					$mut_genotype[1] *= 10**random_normal(1,0,1);
					if ($mut_genotype[1]<=0){
						while ($mut_genotype[1]<=0){
							$mut_genotype[1] = $genotype[1] * 10**random_normal(1,0,1);
						}
					}
					$mutation = "rho";
				}
				if ($mut_rand[2]==1 and $testing_capacitance == 0){
					$mutation = "trait";
					if (rand() > 30/330){
						$trait_mut = 1; $mutation = "alpha";
						my $locus = int($num_trait_loci * rand());
						my $trait = int($num_traits * rand());
						if (rand()<3/300 & $cooption ==1 & $mut_traitdel[$locus] == 0){
						#print "coopt\n";
							$mut_alpha[$locus][$trait] += $mut_beta[$locus][$trait];
							my $beta_var = ($sigmam/$num_trait_loci)**2 / (1 - ((1/$trait_mut_bias -1)*$trait_mut_bias )**2);
							$mut_beta[$locus][$trait] = 0;#random_normal(0, $beta_var);
						}else{
							$mut_alpha[$locus][$trait] += random_normal(1,-$trait_mut_bias * $mut_alpha[$locus][$trait] ,$sigmam/$num_trait_loci);#(1,-$mut_alpha[$locus][$trait]/50,$sigmam/$num_trait_loci)
						}
					}else{
						$trait_mut = 1; $mutation = "beta";
						my $locus = int($num_trait_loci * rand());
						my $trait = int($num_traits * rand());			
						$mut_beta[$locus][$trait] += random_normal(1,-$mut_beta[$locus][$trait]/$beta_mut_bias,$sigmam/$num_trait_loci);#random_normal(1,-$mut_beta[$locus][$trait]/50,$sigmam/$num_trait_loci);
						if ($mut_traitdel[$locus] == 0 && rand() < 0.4 && $ldel_mut_freq > 0){$mut_traitdel[$locus] = 1;#print "del mutation \n";
						}
						elsif ($mut_traitdel[$locus] == 1 && rand() < 0.1 && $ldel_mut_freq > 0){$mut_traitdel[$locus] = 0;#print "benign mutation \n";
						}
					}
					
				}
			}
			my @tempgenotype = ($Ldel[$deme], $rho[$deme]);
			my $WT_fitness = calculate_fitness(\@tempgenotype, \@curr_alpha,\@curr_beta, \@curr_traitdel, \@env_opt, $bumpiness) ;
            if ($WT_fitness <= 0){$WT_fitness = 0.0001;}
            if ($ldel_mutant == 1 && $testing_capacitance == 1){
           		$WT_fitness = 0;
				if($high_rho_bout == 0){
					$Ldel_index = int($Ldel[$deme]/25) + 1;
					if ($Ldel_index > 16){$Ldel_index = 16;}#print $Ldel_index . "\tindex\n";
				}
           		my @new_rho_and_bout_num = get_test_rho(\@tempgenotype, $high_rho_bout, $high_rho_sample_number, \@high_rho_values, $high_rho_bout_number, $using_mode, \@total_high_rho, $Ldel_index);
				#foreach my $a (@new_rho_and_bout_num){print "got rho\t"  . $a . "\t";print "\n";}
           		my $temp_rho =  $new_rho_and_bout_num[0];
				$high_rho_bout = $new_rho_and_bout_num[1];
				$high_rho_bout_number = $new_rho_and_bout_num[2];				
				$high_rho_sample_number = $new_rho_and_bout_num[3];
				@tempgenotype = ($Ldel[$deme], $temp_rho);
				$WT_fitness = calculate_fitness(\@tempgenotype, \@curr_alpha,\@curr_beta, \@curr_traitdel, \@env_opt) ;
				@mut_genotype = ($mut_genotype[0], $temp_rho);
				if(log($temp_rho)/log(10) > -2.5){$num_high_rho_values_used++;
				#	print $Ldel[$deme]. "\t" .  $temp_rho . "\t" . $WT_fitness . "\n";
				}
			}else{
				@mut_genotype = ($mut_genotype[0], $mut_genotype[1]);
			}
			my $mut_fitness = calculate_fitness(\@mut_genotype, \@mut_alpha,\@mut_beta, \@mut_traitdel, \@env_opt, $bumpiness);
			my $pfix = 0;
			my $rel_fitness = 0;
			if ($WT_fitness > 0){
				$rel_fitness = $mut_fitness / $WT_fitness - 1;
				if ($rel_fitness ==0){
					$pfix = 1/$Ne;
				}else{
					$pfix = (1 - exp(-2*$rel_fitness))/(1 - exp(-2*$Ne*$rel_fitness));
				}
			}

			if (rand() < $pfix){
				$Ldel[$deme] = $mut_genotype[0];
				if($ldel_mutant != 1 || $testing_capacitance != 1){
					$rho[$deme] = $mut_genotype[1];
				}
				for (my $i = 0; $i < $num_trait_loci; $i++){
					for (my $j = 0; $j < $num_traits; $j++){					
						$alpha[$deme][$i][$j]= $mut_alpha[$i][$j];
						$beta[$deme][$i][$j] = $mut_beta[$i][$j];
					}
					$traitdel[$deme][$i] = $mut_traitdel[$i];
				}
			}
			my $log_rho = log($rho[$deme])/log(10);
 			my $rho_bin = int(-1 * $log_rho * 20);
			$rho_frequency[$rho_bin]++;
			if ($log_rho > -2.5 && $gen > $burnin){
				$high_rho_bout = 1;
				$high_rho_duration[$rho_record_num]++;#print $high_rho_duration[$rho_record_num] . "\n";
			}elsif($high_rho_bout == 1 && $gen > $burnin){

			}

		}
		
		if ($gen % $environment_duration==0 && $testing_capacitance == 1){
			$Ldel_index = int($Ldel[0]/25) + 1;
			if ($Ldel_index > 16){$Ldel_index = 16;}
			#print $total_bouts[$Ldel_index]/((10000-6000)*2) . "\t";
			#print $high_rho_bout ;
			#print "\t" . int(scalar(@{$high_rho_values[$Ldel_index]})). "\n";
			if($gen > 1 && $high_rho_bout==0 ){#&& rand() < $total_bouts[$Ldel_index]/((55000-6000)*3)
				$high_rho_sample_number=0;
				$high_rho_bout_number = int(rand()*scalar(@{$high_rho_values[$Ldel_index]}));
				$high_rho_bout=1;
				#print $high_rho_bout_number . "\tstarting bout\n";
				if ($Ldel_index > 16){$Ldel_index = 16;}
			}
		}
		if ($gen % $environment_duration==0 && $gen > 1){
			if ($high_rho_bout == 1){				
				my $i = 0;
				foreach (@rho_recorder){
					#if ($_ > -3.5){
					printf RHOOUTPUT $rho_record_num . "\t" . $gen_recorder[$i] . "\t" . $ldel_recorder[$i] . "\t" . $_ . "\t" . $env_recorder[$i] . "\t" . $alpha_recorder[$i] . "\t" . $beta_recorder[$i] . "\n";
				#	}
					$i++;
				}
				$rho_record_num++;
			}
			$rho_record = 0; #$high_rho_bout = 0;
			@gen_recorder =[];@ldel_recorder=[];@rho_recorder=[];@env_recorder=[];@alpha_recorder=[];@beta_recorder=[];
			if($environmental_change == 1){
				my @new_environmental_optimum;
				my $viable = 0;
				my @curr_alpha;
				my @curr_beta;
				my @mut_alpha;
				my @mut_beta;
				my @curr_traitdel;
				
				for (my $i = 0; $i < $num_trait_loci; $i++){
					for (my $j = 0; $j < $num_traits; $j++){
						$curr_alpha[$i][$j] = $alpha[0][$i][$j];
						$curr_beta[$i][$j] = $beta[0][$i][$j];
					}
					$curr_traitdel[$i] = $traitdel[0][$i];
				}
				while ($viable == 0){
					for (my $i = 0; $i < $num_traits; $i++){
						for (my $deme =0; $deme < $num_demes; $deme++){
							$new_environmental_optimum[$deme][$i] = $environmental_optimum[$deme][$i] + random_normal(1,-$env_bias * $environmental_optimum[$deme][$i],$sigma_e);	
						} 
					}	
					my @env_opt;
					for (my $i = 0; $i < $num_traits; $i++){	
						$env_opt[$i] = $new_environmental_optimum[0][$i];
					}
					my @tempgenotype = ($Ldel[0], $rho[0]);
					my $WT_fitness = calculate_fitness(\@tempgenotype, \@curr_alpha,\@curr_beta, \@curr_traitdel, \@env_opt, 0) ;
					if ($WT_fitness > 0){
						@environmental_optimum = @new_environmental_optimum; $viable = 1; 		
					}
				}
			}
		}
		if (rand() < $num_printouts/($generations-$burnin) && $gen > $burnin ){#||$gen ==0
			my $alpha_sum = 0;
			my $beta_sum = 0;
			for(my $i = 0; $i<$num_trait_loci; $i++){
				for(my $j = 0; $j<$num_traits; $j++){
					$alpha_sum += $alpha[0][$i][$j];
					$beta_sum += $beta[0][$i][$j];
				}
			}
			my $mean_Ldel;
			my $mean_rho;
			my $mean_traitdel;
			for (my $deme =0; $deme < $num_demes; $deme++){
				my $sum_traitdel = 0;
				for (my $i; $i<$num_trait_loci; $i++){
					$sum_traitdel+=$traitdel[$deme][$i];
				}
				$mean_traitdel +=$sum_traitdel;
				$mean_Ldel+=$Ldel[$deme]+$sum_traitdel;
				$mean_rho+=$rho[$deme];
			}
			$mean_traitdel/=$num_demes;
			$mean_Ldel/=$num_demes;
			$mean_rho/=$num_demes;
			my $sum_env_opt;
			print $Ne . "\t" . $Ldel_init[0] .  "\t"  . $using_mode.  "\t" . $gen . "\t" . $mean_Ldel . "\t" . log($mean_rho)/log(10) . "\t" . $environmental_optimum[0][0]. "\t" . $alpha_sum . "\t" . $mean_traitdel ."\t" . $beta_sum . "\t" . $num_high_rho_values_used . "\n";
			printf OUTPUT $Ne . "\t" . $Ldel_init[0] .  "\t"  . $using_mode.  "\t" . $gen . "\t" . $mean_Ldel . "\t" . log($mean_rho)/log(10) . "\t" . $environmental_optimum[0][0]. "\t" . $alpha_sum . "\t" . $mean_traitdel ."\t" . $beta_sum . $num_high_rho_values_used ."\n";
			$num_high_rho_values_used =0;
		}
		
		if($gen % ($environment_duration/10000) == 0 && $record_rho == 1 && $gen > $burnin){
			my $alpha_sum = 0;
			my $beta_sum = 0;
			my $sum_trait_del = 0;
			for(my $i = 0; $i<$num_trait_loci; $i++){
				for(my $j = 0; $j<$num_traits; $j++){
					$alpha_sum += $alpha[0][$i][$j];
					$beta_sum += $beta[0][$i][$j];
				}
				$sum_trait_del += $traitdel[$i];
			}
			my $mean_Ldel;
			my $mean_rho = 0;
			my $mean_traitdel;
			for (my $deme =0; $deme < $num_demes; $deme++){
				my $sum_traitdel = 0;
				for (my $i=0; $i<$num_trait_loci; $i++){
					$sum_traitdel+=$traitdel[$deme][$i];
				}
				$mean_traitdel +=$sum_traitdel;
				$mean_Ldel+=$Ldel[$deme]+$sum_traitdel;
				$mean_rho+=$rho[$deme];
			}
			#print $alpha_sum . "\t" . $beta_sum . "\n";
			$mean_traitdel/=$num_demes;
			$mean_Ldel/=$num_demes;
			$mean_rho/=$num_demes;
			$rho_recorder[$rho_record] = log($mean_rho)/ log(10); 
			$alpha_recorder[$rho_record] = $alpha_sum;
			$beta_recorder[$rho_record] = $beta_sum;
			$gen_recorder[$rho_record] = $gen;
			$env_recorder[$rho_record] = $environmental_optimum[0][0];
			$ldel_recorder[$rho_record]  = $mean_Ldel;
			$rho_record++;
		}	
	}	
	my $durationfilename = "high_rho_duration_" . $pop_magnitude . "_" . $Ldel_init[0] . "_" . $using_mode . "_" . $replicate . ".txt";
	open DURATIONOUTPUT, '>', $durationfilename or die "Can't open $filename - $!\n";
	foreach(@high_rho_duration){
		printf DURATIONOUTPUT $_ . "\n";
	}
	close DURATIONOUTPUT;

	my $rho_bin = 0;
	foreach (@rho_frequency){
		$rho_bin -= 0.05;
		printf RHOFREQOUTPUT  $rho_bin . "\t" . $_ . "\n";
	}
	close RHOFREQOUTPUT;
	close RHOOUTPUT;			
}

close OUTPUT;

sub calculate_fitness {#(\@tempgenotype, \@curr_alpha,\@curr_beta, \@curr_traitdel, \@env_opt, 0) ;

	my @genotype = @{$_[0]};
	my @alpha= @{$_[1]};
	my @beta= @{$_[2]};
	my @traitdel= @{$_[3]};
	my @environmental_optimum = @{$_[4]};
	my $sum_trait_del;
	my $curr_bumpiness = $_[5];
	

	foreach(@traitdel){$sum_trait_del+=$_;}
	my $Ldel = $genotype[0]+$sum_trait_del;
	my $rho = $genotype[1];
    my $omega_m = 1 - $gamma * $rho * $Ldel / $num_del_loci;
	my $high_rho_penalty = 1 - ($rho / 0.5)**20;
	my $omega_e = 1;
	my $best_omega_e = 0;
	my @trait_values;
	for (my $j = 0; $j < $num_traits; $j++){
		for (my $i = 0; $i < $num_trait_loci; $i++){
			$trait_values[$j] += $alpha[$i][$j] + $rho * (1-$traitdel[$i]) * $beta[$i][$j];
		}
	}
	for (my $i = 0; $i < ($curr_bumpiness + 1); $i++){
		my $curr_sigma_f = $sigma_f / ($curr_bumpiness + 1);
		my $curr_environmental_optimum = $environmental_optimum[0] + sqrt(-$sigma_f**2 * log(1/(2*$curr_bumpiness+1)));
		#print $curr_bumpiness . "\t" .$rho . "\t" .  $curr_environmental_optimum . "\t";

		for (my $sign = -1; $sign < 2; $sign += 2){
			my $sum_deviations = 0;
			
			for (my $j = 0; $j < $num_traits; $j++){
				#print $trait_values[$j] . "\n";
				$sum_deviations += ($trait_values[$j] - $sign * $curr_environmental_optimum)**2;
			}
			$omega_e = exp((-1*$sum_deviations)/(2*$curr_sigma_f**2));
			if ($best_omega_e < $omega_e or $best_omega_e == 0){
				$best_omega_e = $omega_e;
			}
		}
	}
	$omega_e = $best_omega_e;

    if ($omega_m<0){$omega_m = 0;}
	if ($high_rho_penalty<0){$high_rho_penalty = 0;}

    #my $omega_a = (1 - $mu_a)** $Ldel;
    my $omega_s = 1/(1 - log($rho) * $delta);
	my $fitness = $omega_e * $omega_m * $omega_s *$high_rho_penalty;	
	#print $omega_e . "\t" . $omega_m . "\t" . $omega_s . "\t" . $high_rho_penalty . "\n";
	return $fitness;
	
}

sub change_environment {
	my @curr_opt = @{$_[0]};
	my $env_bias = $_[1];
	my $sigma_e = $_[2];
	my $num_demes = $_[3];
	my $num_traits = $_[4];
	my @new_environmental_optimum;
	for (my $deme = 0; $deme < $num_demes; $deme++){	
		for (my $i = 0; $i < $num_traits; $i++){	
			$new_environmental_optimum[$deme][$i] = $curr_opt[$deme][$i] + random_normal(1,-$env_bias * $curr_opt[$deme][$i],$sigma_e);	
		}
	}
	return @new_environmental_optimum; 
}

sub get_test_rho {
	my @genotype = @{$_[0]};
	my $high_rho_bout = $_[1];
	my $high_rho_sample_number = $_[2];
	my @high_rho_values = @{$_[3]};
	my $high_rho_bout_number = $_[4];	
	my $using_mode = $_[5];
	my @total_high_rho = @{$_[6]};
	my $prob_start_capacitance;
	my $prob_high_ldel;
	my @low_rho_variances=
([1.074749962,0.002984431],
[1.150672396,0.00312084],
[1.204159525,0.003031377],
[1.234705743,0.003141654],
[1.248120966,0.003083763],
[1.248532199,0.003133024],
[1.174932014,0.003124239],
[1.118509457,0.003189425],
[0.80538015,0.002960043],
[0.553665418,0.003124822],
[0.166238768,0.003109595],
[0.112111343,0.003033184],
[0.06483892,0.003112729],
[0.061106424,0.003140483],
[0.048105472,0.003105317],
[0.045181653,0.00307821]);
	my $temp_rho = $genotype[1];
	my $Ldel = $genotype[0];
	my $rho = $genotype[1];
	my $Ldel_index = $_[7];
#	if ($Ldel_index > 16){$Ldel_index = 16;}#print $Ldel_index . "\tindex\n";

	#print $low_rho_variances[$Ldel_index][0] . "\t" . $low_rho_variances[$Ldel_index][0] . "\n";

	if ($using_mode eq "time" && $high_rho_bout==1){
		my @high_rho_bouts = @{$high_rho_values[$Ldel_index]};
		#	print "getting value\t" . scalar(@{$high_rho_bouts[$high_rho_bout_number]}) . "\n";
		if($high_rho_sample_number < scalar(@{$high_rho_bouts[$high_rho_bout_number]})){
			$temp_rho = 10**$high_rho_values[$Ldel_index][$high_rho_bout_number][int($high_rho_sample_number)];
			$high_rho_sample_number+=5000/$environment_duration;
		}else{
			$high_rho_sample_number = 0;
			$high_rho_bout = 0;
		}
	}
	my $prob_high_rho = ($total_high_rho[$Ldel_index] * $environment_duration / 5000)/(15000*(55000-6000)*3);
	if (($using_mode eq "high" || $using_mode eq "both") && rand()<$prob_high_rho){
		$temp_rho = 10**$high_rho_values[$Ldel_index][int(rand()*scalar(@{$high_rho_values[$Ldel_index]}))];#print "high\t$temp_rho\n";
	}elsif($using_mode eq "both" || $using_mode eq "low"){
		$temp_rho = $rho*10**(random_normal(1, 0, sqrt($low_rho_variances[$Ldel_index][0]**2 - $low_rho_variances[$Ldel_index][1]**2)));#print "low\t$temp_rho\n";
	}
	my @returnarray = ($temp_rho, $high_rho_bout, $high_rho_bout_number, $high_rho_sample_number);
	return @returnarray;
}