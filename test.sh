#conda create --name ERICA python=3.6 tensorflow=2.1.0 plotnine=0.6.0
#conda env create --file environment.yml
#source activate ERICA

# four-taxon model
python vcf2MSA.py \
-i test/pop_test.vcf.gz \
-r test/pop_test.fasta \
-o test/pop_test \
-f diplo \
-P1 H_m_aglaope_1,H_m_aglaope_2,H_m_aglaope_3,H_m_aglaope_4 \
-P2 H_m_amaryllis_1,H_m_amaryllis_2,H_m_amaryllis_3,H_m_amaryllis_4 \
-P3 H_t_thelxinoe_1,H_t_thelxinoe_2,H_t_thelxinoe_3,H_t_thelxinoe_4 \
-P4 H_ethilla_1,H_ethilla_2,H_ethilla_3,H_ethilla_4

python ERICAPrediction.py -i test/pop_test_Hmel218003o.txt -o test/pop_test_Hmel218003o_res.txt -p 4 

python ERICAVisualization.py \
-i test/pop_test_Hmel218003o_res.txt \
-o test/pop_test_Hmel218003o_10k \
-p 4 \
-w 10000 \
-r 1:200 \
-c Hmel218003o \
-d 0.5


# five-taxon model 
python vcf2MSA.py \
-i test/five_pop_test.vcf.gz \
-r test/five_pop_test.fasta \
-o test/five_pop_test \
-f haplo \
-P1 GP39,GP77,GP536,GP640,GP761-1 \
-P2 Nipponbare,HP14,HP44,HP48,HP314,HP103,HP45,UR28 \
-P3 W1943,W3095-2,Orufi,W3078-2 \
-P4 W0170,W1698,W1754,W0123-1,Oniva \
-P5 Obart

python ERICAPrediction.py -i test/five_pop_test_10.txt -o test/five_pop_test_10_res.txt -p 5 

python ERICAVisualization.py \
-i test/five_pop_test_10_res.txt \
-o test/five_pop_test_10_10k \
-p 5 \
-w 10000 \
-r 1:200 \
-c Chr10 \
-d 0.5
