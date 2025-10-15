from snakemake.script import snakemake
import logging
import polars as pl
import polars.selectors as cs

# Log file edition
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

correspondence_classification = [

    # N.B.: For WOOD IN THE ROUGH (ROUNDWOOD) at the aggregated level choice is made
    # to allocate 440111.12 to WOOD CHIPS AND PARTICLES, RESIDUES AND RECOVERABLE WOOD PRODUCTS
    # and to allocate 440311.12 to it from WOOD SIMPLY WORKED OR PROCESSED
    # for simplicity and clarity.

    {
        'FAO Product Agg': 'WOOD IN THE ROUGH (ROUNDWOOD)',
        'FAO Code Agg': '01',
        'FAO 1982': ('1',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4403,),
        'HS 2002': (4403,),
        'HS 2007': (4403,),
        'HS 2012': (4403,),
        'HS 2017': (4403,),
        'HS 2022': (4403,)
    },

    {
        'FAO Product Agg': 'WOOD IN THE ROUGH (ROUNDWOOD)',
        'FAO Code Agg': '01',
        'FAO 1982': ('1',),
        'FAO Product': ('Wood fuel (including wood for charcoal)',),
        'FAO Code': ('011',),
        'HS 1996': (440110, 440110),
        'HS 2002': (440110, 440110),
        'HS 2007': (440110, 440110),
        'HS 2012': (440110, 440110),
        'HS 2017': (440111, 440112),
        'HS 2022': (440111, 440112)
    },

    {
        'FAO Product Agg': 'WOOD IN THE ROUGH (ROUNDWOOD)',
        'FAO Code Agg': '01',
        'FAO 1982': ('1',),
        'FAO Product': ('Wood in the rough other than wood fuel',),
        'FAO Code': ('012',),
        'HS 1996': (440320, 440320, 440320, 440320, 440320, 440320, 440341, 440349, 440349, 440391, 440392, 440392, 440399, 440399, 440399, 440399, 440399),
        'HS 2002': (440320, 440320, 440320, 440320, 440320, 440320, 440341, 440349, 440349, 440391, 440392, 440392, 440399, 440399, 440399, 440399, 440399),
        'HS 2007': (440320, 440320, 440320, 440320, 440320, 440320, 440341, 440349, 440349, 440391, 440392, 440392, 440399, 440399, 440399, 440399, 440399),
        'HS 2012': (440320, 440320, 440320, 440320, 440320, 440320,  # -> 440321 | 440322 | 440323 | 440324 | 440325 | 440326 in HS 2017
                    440341, 440349, 440349, 440391, 
                    440392, 440392,  # -> 440393 | 440394 in HS 2017
                    440399, 440399, 440399, 440399, 440399),  # -> 440395 | 440396 | 440397 | 440398 | 440399 in HS 2017
        'HS 2017': (440321, 440322, 440323, 440324, 440325, 440326,  # new in HS 2017
                    440341, 
                    440349, 440349,  # -> 440342 | 440349 in HS 2022
                    440391, 
                    440393, 440394,  # new in HS 2017
                    440395, 440396, 440397, 440398,  # new in HS 2017
                    440399),
        'HS 2022': (440321, 440322, 440323, 440324, 440325, 440326, 440341, 
                    440342, 440349,  # new in HS 2022
                    440391, 440393, 440394, 440395, 440396, 440397, 440398, 440399)
    },

    # N.B.: For 'WOOD SIMPLY WORKED OR PROCESSED', we choose to not include 
    # 4402.XX and 4403.XX at the aggregated level for clarity and simplicity.

    {
        'FAO Product Agg': 'WOOD SIMPLY WORKED OR PROCESSED',
        'FAO Code Agg': '02',
        'FAO 1982': ('4',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4404, 4405),
        'HS 2002': (4404, 4405),
        'HS 2007': (4404, 4405),
        'HS 2012': (4404, 4405),
        'HS 2017': (4404, 4405),
        'HS 2022': (4404, 4405)
    },

    # N.B.: For 'WOOD SIMPLY WORKED OR PROCESSED', product 'Wood fuel simply worked'
    # has not been included as the link between classification is partial and as
    # it duplicate with the product 'Wood fuel (including wood for charcoal)'.

    {
        'FAO Product Agg': 'WOOD SIMPLY WORKED OR PROCESSED',
        'FAO Code Agg': '02',
        'FAO 1982': ('4',),
        'FAO Product': ('Wood charcoal', 'Torrefied wood'),
        'FAO Code': ('021', '022'),
        'HS 1996': (440200,),  # ***
        'HS 2002': (440200,),  # *** - 440200 is subdivided in 440210 and 440290 (440210 of bamboo)
        'HS 2007': (440290,),
        'HS 2012': (440290,),
        'HS 2017': (440290,),
        'HS 2022': (440290,)
    },

    {
        'FAO Product Agg': 'WOOD SIMPLY WORKED OR PROCESSED',
        'FAO Code Agg': '02',
        'FAO 1982': ('4',),
        'FAO Product': ('Roundwood treated with preservatives',),
        'FAO Code': ('023',),
        'HS 1996': (440310, 440310),
        'HS 2002': (440310, 440310),
        'HS 2007': (440310, 440310),
        'HS 2012': (440310, 440310),  # ->  440311 | 440312 in HS 2017
        'HS 2017': (440311, 440312),
        'HS 2022': (440311, 440312)
    },

    {
        'FAO Product Agg': 'WOOD SIMPLY WORKED OR PROCESSED',
        'FAO Code Agg': '02',
        'FAO 1982': ('4',),
        'FAO Product': ('Roughly trimmed wood',),
        'FAO Code': ('024',),
        'HS 1996': (440410, 440420),
        'HS 2002': (440410, 440420),
        'HS 2007': (440410, 440420),
        'HS 2012': (440410, 440420),
        'HS 2017': (440410, 440420),
        'HS 2022': (440410, 440420)
    },

    {
        'FAO Product Agg': 'WOOD SIMPLY WORKED OR PROCESSED',
        'FAO Code Agg': '02',
        'FAO 1982': ('4',),
        'FAO Product': ('Wood wool', 'Wood flour'),
        'FAO Code': ('025', '026'),
        'HS 1996': (440500,),
        'HS 2002': (440500,),
        'HS 2007': (440500,),
        'HS 2012': (440500,),  
        'HS 2017': (440500,),
        'HS 2022': (440500,)
    },

    {
        'FAO Product Agg': 'WOOD CHIPS AND PARTICLES, RESIDUES AND RECOVERABLE WOOD PRODUCTS',
        'FAO Code Agg': '03',
        'FAO 1982': ('2',),
        'FAO Product': ('Wood residues (including wood for agglomerates)',),
        'FAO Code': ('031',),
        'HS 1996': (440130,),  # ***
        'HS 2002': (440130,),  # ***
        'HS 2007': (440130,),  # *** - 440130 subdivided in 440131 | 440139 in HS 2012 (440131 becomes wood pelets)
        'HS 2012': (440139,),  # *** - 440139 subdivided in 440139 | 440140 in HS 2017
        'HS 2017': (440140,),  # -> 440141 | 440149 in HS 2022
        'HS 2022': (440141,)
    },

    {
        'FAO Product Agg': 'WOOD CHIPS AND PARTICLES, RESIDUES AND RECOVERABLE WOOD PRODUCTS',
        'FAO Code Agg': '03',
        'FAO 1982': ('2',),
        'FAO Product': ('Wood residues (including wood for agglomerates)', 'Recoverable wood products'),
        'FAO Code': ('031', '033'),
        'HS 1996': (440130,),  # ***
        'HS 2002': (440130,),  # ***
        'HS 2007': (440130,),  # *** - 440130 subdivided in 440131 | 440139 in HS 2012 (440131 becomes wood pelets)
        'HS 2012': (440139,),  # *** - 440139 subdivided in 440139 | 440140 in HS 2017
        'HS 2017': (440140,),  # -> 440141 | 440149 in HS 2022
        'HS 2022': (440149,)
    },

    {
        'FAO Product Agg': 'WOOD CHIPS AND PARTICLES, RESIDUES AND RECOVERABLE WOOD PRODUCTS',
        'FAO Code Agg': '03',
        'FAO 1982': ('2',),
        'FAO Product': ('Wood chips and particles',),
        'FAO Code': ('032',),
        'HS 1996': (440121, 440122),
        'HS 2002': (440121, 440122),
        'HS 2007': (440121, 440122),
        'HS 2012': (440121, 440122),  
        'HS 2017': (440121, 440122),
        'HS 2022': (440121, 440122)
    },

    {
        'FAO Product Agg': 'WOOD PELLETS AND OTHER AGGLOMERATES',
        'FAO Code Agg': '04',
        'FAO 1982': ('4122 | 4222',),
        'FAO Product': ('Wood pellets',),
        'FAO Code': ('041',),
        'HS 1996': (440130,),  # ***
        'HS 2002': (440130,),  # ***
        'HS 2007': (440130,),  # *** - 440130 subdivided in 440131 and 440139 in HS 2012 (440131 becomes wood pelets)
        'HS 2012': (440131,),  
        'HS 2017': (440131,),
        'HS 2022': (440131,)
    },

    {
        'FAO Product Agg': 'WOOD PELLETS AND OTHER AGGLOMERATES',
        'FAO Code Agg': '04',
        'FAO 1982': ('4122 | 4222',),
        'FAO Product': ('Wood briquettes',),
        'FAO Code': ('042',),
        'HS 1996': (440130,),  # ***
        'HS 2002': (440130,),  # ***
        'HS 2007': (440130,),  # *** - 440130 subdivided in 440131 and 440139 in HS 2012 (440131 becomes wood pelets)
        'HS 2012': (440139,),  # *** - 440139 subdivided in 440139 | 440140 in HS 2017
        'HS 2017': (440139,),  # *** - 440139 subdivided in 440132 | 440139 in HS 2022
        'HS 2022': (440132,)
    },

    {
        'FAO Product Agg': 'WOOD PELLETS AND OTHER AGGLOMERATES',
        'FAO Code Agg': '04',
        'FAO 1982': ('4122 | 4222',),
        'FAO Product': ('Other agglomerates',),
        'FAO Code': ('043',),
        'HS 1996': (440130,),  # ***
        'HS 2002': (440130,),  # ***
        'HS 2007': (440130,),  # *** - 440130 subdivided in 440131 and 440139 in HS 2012 (440131 becomes wood pelets)
        'HS 2012': (440139,),  # *** - 440139 subdivided in 440139 | 440140 in HS 2017
        'HS 2017': (440139,),  # *** - 440139 subdivided in 440132 | 440139 in HS 2022
        'HS 2022': (440139,)
    },

    # N.B.: In the case of sawnwood, one must differenciate the aggregated product
    # that does not specificy any FAO product (4406, 4407), and the FAO products
    # that bring details at the HS-6 digit codes.
    
    {
        'FAO Product Agg': 'SAWNWOOD',
        'FAO Code Agg': '05',
        'FAO 1982': ('5',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4406, 4407),  # ***
        'HS 2002': (4406, 4407),
        'HS 2007': (4406, 4407),
        'HS 2012': (4406, 4407),
        'HS 2017': (4406, 4407),
        'HS 2022': (4406, 4407)
    },

    # Between HS 2017 and HS 2022, products {440711, 440712, 440719} are subdivided 
    # to create new subheadings {440713, 440714} alongside them. We choose to as 
    # a filiation to {440713, 440714} the product 440719 'Other' for simplicity. 
    
    {
        'FAO Product Agg': 'SAWNWOOD',
        'FAO Code Agg': '05',
        'FAO 1982': ('5',),
        'FAO Product': ('Coniferous, sawnwood and sleepers',),
        'FAO Code': ('051',),
        'HS 1996': (440610, 440690,  # ***
                    440710, 440710, 440710, 440710, 440710),  
        'HS 2002': (440610, 440690,  # ***
                    440710, 440710, 440710, 440710, 440710),  
        'HS 2007': (440610, 440690,  # ***
                    440710, 440710, 440710, 440710, 440710),  
        'HS 2012': (440610, 440690,  # *** - 440610 | 440690 subdivided in {440611 | 440612} | {440691 | 440692} in HS 2012
                    440710, 440710, 440710, 440710, 440710),  # 440710 subdivided in 440711 | 440712 | 440719 in HS 2017
        'HS 2017': (440611, 440691, 
                    440711, 440712, 440719, 440719, 440719),  # 440711 | 440712 | 440719 subdivided in 440711 | 440712 | 440713 | 440714 | 440719 in HS 2022
        'HS 2022': (440611, 440691, 440711, 440712, 440713, 440714, 440719)
    },
    
    {
        'FAO Product Agg': 'SAWNWOOD',
        'FAO Code Agg': '05',
        'FAO 1982': ('5',),
        'FAO Product': ('Tropical non-coniferous, sawnwood and sleepers',),
        'FAO Code': ('052',),
        'HS 1996': (440610, 440690,  # ***
                    440724, 440724, 440729, 440725, 440726, 440729, 440729, 440729),  # *** - ex440799 and 449729 partially merged in 440729 in HS 2002
        'HS 2002': (440610, 440690,  # ***
                    440724, 440724, 440729, 440725, 440726, # 440724 subdivided in 440721 | 440722 in HS 2007
                    440729, 440729, 440729),  # 440729 subdivided in 440727 | 440728 | 440729 in HS 2007
        'HS 2007': (440610, 440690,  # ***
                    440721, 440722, 440729, 440725, 440726, 440727, 440728, 440729),  
        'HS 2012': (440610, 440690,  # *** - 440610 | 440690 subdivided in {440611 | 440612} | {440691 | 440692} in HS 2012
                    440721, 440722, 440729, 440725, 440726, 440727, 440728, 440729),  # *** - deletion of Subheading Note 2 to Chapter 44 results in the expansion of the scope of the expression "tropical wood".
        'HS 2017': (440612, 440692,  # ***
                    440721, 440722, 440729, 440725, 440726, 440727, 440728, 440729),  # 440729 subdivided in 440723 | 440729 in HS 2022
        'HS 2022': (440612, 440692,  # ***
                    440721, 440722, 440723, 440725, 440726, 440727, 440728, 440729)
    },

    {
        'FAO Product Agg': 'SAWNWOOD',
        'FAO Code Agg': '05',
        'FAO 1982': ('5',),
        'FAO Product': ('Other non-coniferous, sawnwood and sleepers',),
        'FAO Code': ('053',),
        'HS 1996': (440612, 440692,  # ***
                    440791, 440792, 
                    440799, 440799, 440799, 440799, 440799, 440799),  # *** - ex440799 - the types of tropical wood of subheadings 4403.41 to 4403.49, 4407.24 to 4407.29, 4408.31 to 4408.39 and 4412.13 to 4412.99 have been added to Subheading Note 1 to Chapter 44.
        'HS 2002': (440612, 440692,  # ***
                    440791, 440792, 
                    440799, 440799, 440799, # 440799 subdivided in 440793 | 440794 | 440795 in HS 2007
                    440799, 440799, 440799),  # ***  
        'HS 2007': (440612, 440692,  # ***
                    440791, 440792, 440793, 440794, 440795, 
                    440799, 440799,  # 440799 subdivided in 440796 | 440797 | 440799 in HS 2012
                    440799),  # ***
        'HS 2012': (440612, 440692,  # ***
                    440791, 440792, 440793, 440794, 440795, 440796, 440797, 
                    440799),  # *** - from ex440799 - deletion of Subheading Note 2 to Chapter 44 results in the expansion of the scope of the expression "tropical wood".
        'HS 2017': (440612, 440692,  # ***
                    440791, 440792, 440793, 440794, 440795, 440796, 440797, 440799),  
        'HS 2022': (440612, 440692,  # ***
                    440791, 440792, 440793, 440794, 440795, 440796, 440797, 440799)
    },

    {
        'FAO Product Agg': 'VENEER SHEETS',
        'FAO Code Agg': '06',
        'FAO 1982': ('512',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4408,),
        'HS 2002': (4408,),
        'HS 2007': (4408,),
        'HS 2012': (4408,),
        'HS 2017': (4408,),
        'HS 2022': (4408,)
    },

    {
        'FAO Product Agg': 'VENEER SHEETS',
        'FAO Code Agg': '06',
        'FAO 1982': ('512',),
        'FAO Product': ('Decorative veneer sheets',),
        'FAO Code': ('061',),
        'HS 1996': (440810, 440831, 440839, 440890),  # ***
        'HS 2002': (440810, 440831, 440839, 440890),  # ***  
        'HS 2007': (440810, 440831, 440839, 440890),  # *** 
        'HS 2012': (440810, 440831, 440839, 440890),  # *** 
        'HS 2017': (440810, 440831, 440839, 440890),  # ***
        'HS 2022': (440810, 440831, 440839, 440890)  # ***
    },

    {
        'FAO Product Agg': 'VENEER SHEETS',
        'FAO Code Agg': '06',
        'FAO 1982': ('512',),
        'FAO Product': ('Non-decorative veneer sheets',),
        'FAO Code': ('062',),
        'HS 1996': (440810, 440831, 440839, 440890),  # ***
        'HS 2002': (440810, 440831, 440839, 440890),  # ***  
        'HS 2007': (440810, 440831, 440839, 440890),  # *** 
        'HS 2012': (440810, 440831, 440839, 440890),  # *** 
        'HS 2017': (440810, 440831, 440839, 440890),  # ***
        'HS 2022': (440810, 440831, 440839, 440890)  # ***
    },

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4410, 4411, 4413, 6808),
        'HS 2002': (4410, 4411, 4413, 6808),
        'HS 2007': (4410, 4411, 4413, 6808),
        'HS 2012': (4410, 4411, 4413, 6808),
        'HS 2017': (4410, 4411, 4413, 6808),
        'HS 2022': (4410, 4411, 4413, 6808)
    },

    # N.B.1: In the case of wood-based panels, code 441899 is not taken into account 
    # because of partial link and because it is not part of the wood products of 
    # this category (it only appears in the HS codes of the aggregated products).

    # N.B.2: Concerning plywood, complete restructuration of HS Heading was operated
    # between HS 2002 and HS 2007. HS-6 digit product filiation cannot be maintained
    # for HS 1996 and HS 2002.
    # Veneer and plywood cannot be differentiated in HS 1996.

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Plywood',),
        'FAO Code': ('071',),
        'HS 1996': (4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412, 4412),  # *** - restructuration of 4403, 4408, 4412
        'HS 2002': (4412, 4412, 4412, 4412,  # *** -  441213 | 441214 | 441219 become 441210 (bamboo) | 441231 | 441232 | 441239 in HS 2007
                    4412, 4412, 4412,  # *** - 441299 dispateched between 441210.94.99 in HS 2007
                    4412, 4412, 4412,  # *** - ex441222.29.92.99 partially merged into 441294 in HS 2007
                    4412, 4412, 4412),  # ***
        'HS 2007': (441231, 441232, 441232, 441239, 441299, 441299, 441299, 441294, 441294, 441294, 441299, 441299, 441299),  
        'HS 2012': (441231, 
                    441232, 441232,  # ex441232 subdivided in 441233 | 441234 in HS 2017 (dispatched in 441231.33.34)
                    441239, 
                    441299, 441299, 441299, 441294, 441294, 441294, 441299, 441299, 441299),
        'HS 2017': (441231, 441233, 441234, 441239, 
                    441299, 441299, 441299,  # ex441299 subdivided in 441241 | 441242 | 441249 in HS 2022
                    441294, 441294, 441294,  # ex441294 subdivided in 441251 | 441252 | 441259 in HS 2022
                    441299, 441299, 441299),  # ex441299 subdivided in 441291 | 441292 | 441299 in HS 2022
        'HS 2022': (441231, 441233, 441234, 441239, 441241, 441242, 441249, 441251, 441252, 441259, 441291, 441292, 441299)
    },

    # N.B.: Concerning particle board, complete restructuration of HS Heading was opereted
    # between HS 2002 and HS 2007. HS-6 digit product filiation cannot be maintained
    # for HS 1996, HS 2002, and HS 2007.

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Particle board',),
        'FAO Code': ('072',),
        'HS 1996': (4410, 4410, 4410),  # *** - 441011.19 subdivided to create 441021.29.32.33.39 in HS 2022
        'HS 2002': (4410,  # *** - 441031.32.33.39 partially merged into 441011 in HS 2007
                    4410,  # *** - 441021.29.32.33.39 partially merged into 441019 in HS 2007
                    4410),  
        'HS 2007': (441011, 441019, 441090),  
        'HS 2012': (441011, 441019, 441090),
        'HS 2017': (441011, 441019, 441090),  
        'HS 2022': (441011, 441019, 441090)
    },

    # N.B.: Concerning oriented strand board, complete restructuration of HS Heading 
    # was opereted between HS 2002 and HS 2007. HS-6 digit product filiation cannot 
    # be maintained for HS 1996, HS 2002, and HS 2007.

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Oriented strand board (OSB)',),
        'FAO Code': ('073',),
        'HS 1996': (4410,),  # *** - 441011.19 partially merged to create 441021.29 in HS 2002  
        'HS 2002': (4410,),  # *** - 441021.29 partially merged to create 441012 in HS 2007  
        'HS 2007': (441012,),  
        'HS 2012': (441012,),
        'HS 2017': (441012,),  
        'HS 2022': (441012,)
    },

    # N.B.: Concerning fibreboard, complete restructuration of HS Heading 
    # was opereted between HS 2002 and HS 2007. HS-6 digit product filiation cannot 
    # be maintained for HS 1996, HS 2002, and HS 2007.

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Fibreboard',),
        'FAO Code': ('074',),
        'HS 1996': (4411, 4411, 4411, 4411, 4411, 4411),  # *** 
        'HS 2002': (4411, 4411, 4411, 4411, 4411, 4411),  # *** - complete restructuration of heading 4411 in HS 2007
        'HS 2007': (441112, 441113, 441114, 441192, 441193, 441194),  
        'HS 2012': (441112, 441113, 441114, 441192, 441193, 441194),
        'HS 2017': (441112, 441113, 441114, 441192, 441193, 441194),  
        'HS 2022': (441112, 441113, 441114, 441192, 441193, 441194)
    },

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Densified wood',),
        'FAO Code': ('075',),
        'HS 1996': (441300,),  
        'HS 2002': (441300,),  
        'HS 2007': (441300,),  
        'HS 2012': (441300,),
        'HS 2017': (441300,),  
        'HS 2022': (441300,)
    },

    # N.B.: Concerning comnination board, complete restructuration of HS Heading 
    # was opereted between HS 2002 and HS 2007. HS-6 digit product filiation cannot 
    # be maintained for HS 1996 and HS 2002.

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Combination board',),
        'FAO Code': ('076',),
        'HS 1996': (4410,),  # *** - 441011.19 partially subdivided in 441021.29.31.32.33.39 in HS 2002
        'HS 2002': (4410,),  # *** - 441021.29.31.32.33.39 partially merged to create 441019 in HS 2002
        'HS 2007': (441019,),  # ***
        'HS 2012': (441019,),  # ***
        'HS 2017': (441019,),  # ***
        'HS 2022': (441019,)  # ***
    },

    {
        'FAO Product Agg': 'WOOD-BASED PANELS',
        'FAO Code Agg': '07',
        'FAO 1982': ('6',),
        'FAO Product': ('Other panels based on wood or other 65 n.a. ligno-cellulosic materials',),
        'FAO Code': ('077',),
        'HS 1996': (680800,),  
        'HS 2002': (680800,),  
        'HS 2007': (680800,),  
        'HS 2012': (680800,),
        'HS 2017': (680800,),  
        'HS 2022': (680800,)
    },

    {
        'FAO Product Agg': 'WOOD PULP',
        'FAO Code Agg': '08',
        'FAO 1982': ('7',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4701, 4702, 4703, 4704, 4705),  
        'HS 2002': (4701, 4702, 4703, 4704, 4705),
        'HS 2007': (4701, 4702, 4703, 4704, 4705),
        'HS 2012': (4701, 4702, 4703, 4704, 4705),
        'HS 2017': (4701, 4702, 4703, 4704, 4705),
        'HS 2022': (4701, 4702, 4703, 4704, 4705)
    },

    {
        'FAO Product Agg': 'WOOD PULP',
        'FAO Code Agg': '08',
        'FAO 1982': ('7',),
        'FAO Product': ('Mechanical wood pulp',),
        'FAO Code': ('081',),
        'HS 1996': (470100,),  
        'HS 2002': (470100,),
        'HS 2007': (470100,),
        'HS 2012': (470100,),
        'HS 2017': (470100,),
        'HS 2022': (470100,)
    },

    {
        'FAO Product Agg': 'WOOD PULP',
        'FAO Code Agg': '08',
        'FAO 1982': ('7',),
        'FAO Product': ('Semi-chemical wood pulp',),
        'FAO Code': ('082',),
        'HS 1996': (470500,),  
        'HS 2002': (470500,),
        'HS 2007': (470500,),
        'HS 2012': (470500,),
        'HS 2017': (470500,),
        'HS 2022': (470500,)
    },

    {
        'FAO Product Agg': 'WOOD PULP',
        'FAO Code Agg': '08',
        'FAO 1982': ('7',),
        'FAO Product': ('Chemical wood pulp',),
        'FAO Code': ('083',),
        'HS 1996': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429),  
        'HS 2002': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429),
        'HS 2007': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429),
        'HS 2012': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429),
        'HS 2017': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429),
        'HS 2022': (470311, 470319, 470321, 470329, 470411, 470419, 470421, 470429)
    },

    {
        'FAO Product Agg': 'WOOD PULP',
        'FAO Code Agg': '08',
        'FAO 1982': ('7',),
        'FAO Product': ('Dissolving wood pulp',),
        'FAO Code': ('084',),
        'HS 1996': (470200,),  
        'HS 2002': (470200,),
        'HS 2007': (470200,),
        'HS 2012': (470200,),
        'HS 2017': (470200,),
        'HS 2022': (470200,)
    },

    # N.B.1: For OTHER PULP, code 4702 is not taken into account because of partial
    # link and to not create duplicates with WOOD PULP codes.

    # N.B.2: For pulp from fibers other than wood, complete restructuration of heading to 
    # identify bamboo between HS 2002 and 2007. HS-6 digit product filiation cannot 
    # be maintained for HS 2002 and HS 2007.

    {
        'FAO Product Agg': 'OTHER PULP',
        'FAO Code Agg': '09',
        'FAO 1982': ('77',),
        'FAO Product': ('Pulp from fibres other than wood',),
        'FAO Code': ('091',),
        'HS 1996': (4706, 4706, 4706, 4706, 4706),  # ***
        'HS 2002': (4706, 4706, 4706, 4706, 4706),  # *** - 470691.92.93 subdivided in 470630.91.92.93 in HS 2007
        'HS 2007': (470610, 470630, 470691, 470692, 470693),
        'HS 2012': (470610, 470630, 470691, 470692, 470693),
        'HS 2017': (470610, 470630, 470691, 470692, 470693),
        'HS 2022': (470610, 470630, 470691, 470692, 470693)
    },

    {
        'FAO Product Agg': 'OTHER PULP',
        'FAO Code Agg': '09',
        'FAO 1982': ('77',),
        'FAO Product': ('Recovered fibre pulp',),
        'FAO Code': ('092',),
        'HS 1996': (470620,),
        'HS 2002': (470620,),
        'HS 2007': (470620,),
        'HS 2012': (470620,),
        'HS 2017': (470620,),
        'HS 2022': (470620,)
    },

    {
        'FAO Product Agg': 'RECOVERED PAPER',
        'FAO Code Agg': '10',
        'FAO 1982': ('9',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4707,),  # ***
        'HS 2002': (4707,),
        'HS 2007': (4707,),
        'HS 2012': (4707,),
        'HS 2017': (4707,),
        'HS 2022': (4707,)
    },

    {
        'FAO Product Agg': 'RECOVERED PAPER',
        'FAO Code Agg': '10',
        'FAO 1982': ('9',),
        'FAO Product': ('Unbleached kraft paper or paperboard or corrugated paper or paperboard',),
        'FAO Code': ('101',),
        'HS 1996': (470710,),
        'HS 2002': (470710,),
        'HS 2007': (470710,),
        'HS 2012': (470710,),
        'HS 2017': (470710,),
        'HS 2022': (470710,)
    },

    {
        'FAO Product Agg': 'RECOVERED PAPER',
        'FAO Code Agg': '10',
        'FAO 1982': ('9',),
        'FAO Product': ('Other paper or paperboard made mainly of bleached chemical pulp, not coloured in the mass',),
        'FAO Code': ('102',),
        'HS 1996': (470720,),
        'HS 2002': (470720,),
        'HS 2007': (470720,),
        'HS 2012': (470720,),
        'HS 2017': (470720,),
        'HS 2022': (470720,)
    },

    {
        'FAO Product Agg': 'RECOVERED PAPER',
        'FAO Code Agg': '10',
        'FAO 1982': ('9',),
        'FAO Product': ('Paper or paperboard made mainly of mechanical pulp (newspapers, journals and similar)',),
        'FAO Code': ('103',),
        'HS 1996': (470730,),
        'HS 2002': (470730,),
        'HS 2007': (470730,),
        'HS 2012': (470730,),
        'HS 2017': (470730,),
        'HS 2022': (470730,)
    },

    {
        'FAO Product Agg': 'RECOVERED PAPER',
        'FAO Code Agg': '10',
        'FAO 1982': ('9',),
        'FAO Product': ('Other, including unsorted waste and scrap',),
        'FAO Code': ('109',),
        'HS 1996': (470790,),
        'HS 2002': (470790,),
        'HS 2007': (470790,),
        'HS 2012': (470790,),
        'HS 2017': (470790,),
        'HS 2022': (470790,)
    },
    
    {
        'FAO Product Agg': 'PAPER AND PAPERBOARD',
        'FAO Code Agg': '11',
        'FAO 1982': ('8',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813),  # ***
        'HS 2002': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813),
        'HS 2007': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813),
        'HS 2012': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813),
        'HS 2017': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813),
        'HS 2022': (4801, 4802, 4803, 4804, 4805, 4806, 4808, 4809, 4810, 4812, 4813)
    },

    # N.B.1: For Graphic papers, we do not include subheading 480230 present in 
    # HS 2002 and 1996 but deleted in HS 2007 due to low volume of trade.
    # HS-6 digit product filiation cannot be maintained for HS 1996 and HS 2002.

    {
        'FAO Product Agg': 'PAPER AND PAPERBOARD',
        'FAO Code Agg': '11',
        'FAO 1982': ('8',),
        'FAO Product': ('Graphic papers',),
        'FAO Code': ('111',),
        'HS 1996': (4801,  # *** - the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    4802, 
                    4802,  # *** - 480220 | 480560.70.80 | 482359 restructured into 480220 in HS 2002
                    4802,  # *** - 480251 | 480521.22.23.29.60 |  482359 restructured into 480254 in HS 2002
                    4802, 4802, 4802,  # *** - 480252 | 480521.22.23.29.60 |  482359 restructured into 480255.56.57 in HS 2002
                    4802,  # *** - 480253 | 480521.22.23.29.70 |  482359 restructured into 480258 in HS 2002
                    4802, 4802, 4802,  # *** - 4801 | 480251.52.53.60 | 480521.22.23.29.60.70.80 | 482359 into 480261.62.69 in HS 2002
                    4809, 
                    4810, 4810, 4810,  # *** - 481011.12 | 482340.51.59 restructured into 481013.14.19 in HS 2002
                    4810,  # *** - 481021 | 482340.51.59 restructured into 481022 in HS 2002
                    4810),  # *** - 481011.12.29.39 | 482340.51.59 restructured into 481029 in HS 2002
        'HS 2002': (480100,  # ***
                    480210, 480220, 480254, 480255, 480256, 480257, 480258, 
                    480261, 480262,  # ***
                    480269, 4809, 481013, 481014, 481019, 481022, 481029),
        'HS 2007': (480100,  # ***
                    480210, 480220, 480254, 480255, 480256, 480257, 480258, 
                    480261, 480262,  # ***
                    480269, 4809, 481013, 481014, 481019, 481022, 481029),
        'HS 2012': (480100,  # *** - 480261.62 partially merged with 4801
                    480210, 480220, 480254, 480255, 480256, 480257, 480258, 
                    480261, 480262,  # *** - 480261.62 partially merged with 4801
                    480269, 4809, 481013, 481014, 481019, 481022, 481029),
        'HS 2017': (480100, 480210, 480220, 480254, 480255, 480256, 480257, 480258, 480261, 480262, 480269, 4809, 481013, 481014, 481019, 481022, 481029),
        'HS 2022': (480100, 480210, 480220, 480254, 480255, 480256, 480257, 480258, 480261, 480262, 480269, 4809, 481013, 481014, 481019, 481022, 481029)
    },

    {
        'FAO Product Agg': 'PAPER AND PAPERBOARD',
        'FAO Code Agg': '11',
        'FAO 1982': ('8',),
        'FAO Product': ('Sanitary and household papers',),
        'FAO Code': ('112',),
        'HS 1996': (480300,),
        'HS 2002': (480300,),
        'HS 2007': (480300,),
        'HS 2012': (480300,),
        'HS 2017': (480300,),
        'HS 2022': (480300,)
    },
    {
        'FAO Product Agg': 'PAPER AND PAPERBOARD',
        'FAO Code Agg': '11',
        'FAO 1982': ('8',),
        'FAO Product': ('Packaging materials',),
        'FAO Code': ('113',),
        'HS 1996': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459,  # *** - restructuration since the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    4805, 4805, 4805, 4805, 4805, 4805, 4805, 4805, 4805,  # *** - restructuration since the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    480610, 480620, 480640,  # *** - restructuration since the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    480810, 480820, 480890,  # *** - restructuration since the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    4810, 4810, 4810, 4810, 4810,  # *** - restructuration since dimensional criteria have been deleted for heading 48.10.
                    481131, 481139),  # *** - restructuration amendments regarding the size of paper, paperboard, cellulose wadding and webs of cellulose fibres in strips or rolls
        'HS 2002': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459, 
                    480511, 480512, 480519, 480524, 480525, 480530, 480591, 480592, 480593,
                    480610, 480620, 480640,
                    480810, 480820, 480890,
                    481031, 481032, 481039, 481092, 481099,
                    481151, 481159),  # *** - deletion of heading 48.15 because of the low volume of trade -> included in 481151.59
        'HS 2007': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459, 
                    480511, 480512, 480519, 480524, 480525, 480530, 480591, 480592, 480593,
                    480610, 480620, 480640,
                    480810, 480820, 480890,  # subheadings 480820.30 merged into 480840 because low volume of trade in HS 2012
                    481031, 481032, 481039, 481092, 481099,
                    481151, 481159),
        'HS 2012': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459, 
                    480511, 480512, 480519, 480524, 480525, 480530, 480591, 480592, 480593,
                    480610, 480620, 480640,
                    480810, 480840, 480890,
                    481031, 481032, 481039, 481092, 481099,
                    481151, 481159),
        'HS 2017': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459, 
                    480511, 480512, 480519, 480524, 480525, 480530, 480591, 480592, 480593, 
                    480610, 480620, 480640, 
                    480810, 480840, 480890,
                    481031, 481032, 481039, 481092, 481099,
                    481151, 481159),
        'HS 2022': (480411, 480419, 480421, 480429, 480431, 480439, 480442, 480449, 480451, 480452, 480459, 
                    480511, 480512, 480519, 480524, 480525, 480530, 480591, 480592, 480593, 
                    480610, 480620, 480640, 
                    480810, 480840, 480890,
                    481031, 481032, 481039, 481092, 481099,
                    481151, 481159)
    },

    {
        'FAO Product Agg': 'PAPER AND PAPERBOARD',
        'FAO Code Agg': '11',
        'FAO 1982': ('8',),
        'FAO Product': ('Other paper and paperboard',),
        'FAO Code': ('114',),
        'HS 1996': (480240,  # *** - 480240 | 480560.70 | 482359 partially restructured into 480240 since the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    480441, 480540, 480550, 480630,  # *** - the scope of paper and paperboard of a kind used for writing, printing or other graphic purposes has been expanded
                    481200, 481310, 481320),
        'HS 2002': (480240, 480441, 480540, 480550, 480630, 481200, 481310, 481320),
        'HS 2007': (480240, 480441, 480540, 480550, 480630, 481200, 481310, 481320),
        'HS 2012': (480240, 480441, 480540, 480550, 480630, 481200, 481310, 481320),
        'HS 2017': (480240, 480441, 480540, 480550, 480630, 481200, 481310, 481320),
        'HS 2022': (480240, 480441, 480540, 480550, 480630, 481200, 481310, 481320)
    },

    {
        'FAO Product Agg': 'CORK',
        'FAO Code Agg': '12',
        'FAO 1982': ('10',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4501, 4502, 4503, 4504),  # ***
        'HS 2002': (4501, 4502, 4503, 4504),
        'HS 2007': (4501, 4502, 4503, 4504),
        'HS 2012': (4501, 4502, 4503, 4504),
        'HS 2017': (4501, 4502, 4503, 4504),
        'HS 2022': (4501, 4502, 4503, 4504)
    },

    {
        'FAO Product Agg': 'CORK',
        'FAO Code Agg': '12',
        'FAO 1982': ('10',),
        'FAO Product': ('Natural cork',),
        'FAO Code': ('121',),
        'HS 1996': (450110, 450190, 450200, 450310, 450390),
        'HS 2002': (450110, 450190, 450200, 450310, 450390),
        'HS 2007': (450110, 450190, 450200, 450310, 450390),
        'HS 2012': (450110, 450190, 450200, 450310, 450390),
        'HS 2017': (450110, 450190, 450200, 450310, 450390),
        'HS 2022': (450110, 450190, 450200, 450310, 450390)
    },

    {
        'FAO Product Agg': 'CORK',
        'FAO Code Agg': '12',
        'FAO 1982': ('10',),
        'FAO Product': ('Agglomerated cork and articles of agglomerated cork',),
        'FAO Code': ('122',),
        'HS 1996': (450410, 450490),
        'HS 2002': (450410, 450490),
        'HS 2007': (450410, 450490),
        'HS 2012': (450410, 450490),
        'HS 2017': (450410, 450490),
        'HS 2022': (450410, 450490)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('13',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (440910, 440920, 440920,  # ***
                    4414, 4415, 4416, 4417,
                    441810, 441810,  # preserve length 
                    441820, 441820,  # preserve length
                    441890,
                    441840, 441850,
                    441830, 441830, 441830,
                    441890, 441890, 441890, 441890, 441890,  # *** - preserve length
                    4419, 4419,
                    4420, 442110,
                    442190,
                    442190,  # preserve length
                    940130, 940140,  # ***
                    940161, 940169,
                    940190,  #***
                    940330, 940340, 940350, 940360,
                    940390,  # ***
                    9406),  # ***
        'HS 2002': (440910, 440920, 440920,  # *** - 440920 divided into 440921.29 in HS 2007 (440921 of bamboo)
                    4414, 4415, 4416, 4417,
                    441810, 441810,  # preserve length 
                    441820, 441820,  # preserve length
                    441890,  # -> 441860 in HS 2007
                    441840, 441850,
                    441830, 441830, 441830,  # 441830.90 -> 4418.7X in HS 2007 
                    441890, 441890, 441890, 441890, 441890,  # *** - preserve length
                    4419, 4419,
                    4420, 442110,
                    442190,
                    442190,  # preserve length
                    940130, 940140,  # ***
                    940161, 940169,
                    940190,  #***
                    940330, 940340, 940350, 940360,
                    940390,  # ***
                    9406),  # ***
        'HS 2007': (440910, 440929, 440929,  # preserve length
                    4414, 4415, 4416, 4417,
                    441810, 441810,  # preserve length 
                    441820, 441820,  # preserve length
                    441860,
                    441840, 441850, 
                    441871, 441872, 441879, 
                    441890, 441890, 441890, 441890,  # preserve length
                    441890,  # ***
                    4419, 4419,  # preserve length
                    4420, 442110,
                    442190,
                    442190,  # preserve length
                    940130, 940140,  # ***
                    940161, 940169,
                    940190,  #***
                    940330, 940340, 940350, 940360,
                    940390,  # ***
                    9406),  # ***
        'HS 2012': (440910, 440929, 440929,  # -> 440929 subdivided in 440922 | 440929 in HS 2017
                    4414, 4415, 4416, 4417,
                    441810, 441810,  # preserve length 
                    441820, 441820,  # preserve length
                    441860,
                    441840, 441850, 
                    441871, 441872, 441879,  # -> 441874 | 441875 | 441879 in HS 2017 (mixed with new 441873)
                    441890, 441890, 441890, 441890,  # *** - -> 441890 is subdivided in 441891 | 441899 in HS 2017
                    441890,  # preserve length
                    4419, 4419,  # -> 4419 is subdivided in 441911 | 441912 | 441919 | 441990 in HS 2017 (441990 of wood, rest of bamboo)
                    4420, 442110, 
                    442190,  # -> 442190 subdivided in 442191 | 442199 in HS 2017 (442199 of other, 442191 of bamboo)
                    442190,  # preserve length
                    940130, 940140,  # ***
                    940161, 940169,
                    940190,  #***
                    940330, 940340, 940350, 940360,
                    940390,  # ***
                    9406),  # *** - -> 940600 subdivided in 940610 and 940690 in HS 2017 (940610 of wood)
        'HS 2017': (440910, 440922, 440929,
                    4414, 4415, 4416, 4417,
                    441810, 441810,  # -> 441811 | 441819 in HS 2022
                    441820, 441820,  # -> 441821 | 441829 in HS 2022
                    441860,  # -> 441830 in HS 2022
                    441840, 441850, 
                    441874, 441875, 441879,
                    441899, 441899, 441899, 441899,  # -> 441881 | 441882 | 441883 | 441889 in HS 2022 (mixed with 441860 and 441891)
                    441899,  # ***
                    441990, 441990,  # -> 441920 | 441990 in HS 2022
                    4420, 442110, 
                    442199,
                    442199,  # -> 442120 in HS 2022 (mixed with 442191)
                    940130, 940140,  # *** - -> 940131 | 940141 in HS 2022
                    940161, 940169,
                    940190,  # *** - 940191 in HS 2022
                    940330, 940340, 940350, 940360,
                    940390,  # *** - 940391 in HS 2022
                    940610),
        'HS 2022': (440910, 440922, 440929,
                    4414, 4415, 4416, 4417,
                    441811, 441819,
                    441821, 441829, 
                    441830, 
                    441840, 441850, 
                    441874, 441875, 441879,
                    441881, 441882, 441883, 441889,
                    441899,  # ***
                    441920, 441990,
                    4420, 442110, 442199,
                    442120, 
                    940131, 940141, 
                    940161, 940169, 
                    940191,
                    940330, 940340, 940350, 940360, 
                    940391, 
                    940610)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Further processed sawnwood',),
        'FAO Code': ('131',),
        'HS 1996': (440910, 
                    440920, 440920),  # ***
        'HS 2002': (440910, 
                    440920, 440920),  # *** - 440920 subdivided in 440921.29 in HS 2007 (440921 of bamboo)
        'HS 2007': (440910, 440929, 440929),
        'HS 2012': (440910, 
                    440929, 440929),  # 440929 subdivided in 440922.29 in HS 2017
        'HS 2017': (440910, 440922, 440929),
        'HS 2022': (440910, 440922, 440929)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Wooden wrapping and packaging material',),
        'FAO Code': ('132',),
        'HS 1996': (441510, 441520, 441600),
        'HS 2002': (441510, 441520, 441600),
        'HS 2007': (441510, 441520, 441600),
        'HS 2012': (441510, 441520, 441600),
        'HS 2017': (441510, 441520, 441600),
        'HS 2022': (441510, 441520, 441600)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Wood products for domestic/decorative use',),
        'FAO Code': ('133',),
        'HS 1996': (441400, 441400, 
                    441900, 441900,  # ***
                    442010, 442010, 442090),
        'HS 2002': (441400, 441400, 
                    441900, 441900,  # ***
                    442010, 442010, 442090),
        'HS 2007': (441400, 441400, 
                    441900, 441900,  # ***
                    442010, 442010, 442090),
        'HS 2012': (441400, 441400, 
                    441900, 441900,  # *** - 4419 has been subdivided to provide separately for certain articles of bamboo
                    442010, 442010, 442090),
        'HS 2017': (441400, 441400,  # 441400 subdivided into 441410.90 in HS 2022
                    441990, 441990,  # 441990 subdivided into 441920.90 in HS 2022
                    442010, 442010,  # 442010 subdivided into 442011.19 in HS 2022
                    442090),
        'HS 2022': (441410, 441490, 441920, 441990, 442011, 442019, 442090)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Other manufactured wood products',),
        'FAO Code': ('134',),
        'HS 1996': (441700, 442110, 442110,
                    442190),  # ***
        'HS 2002': (441700, 442110, 442110,
                    442190),  # ***
        'HS 2007': (441700, 442110, 442110,
                    442190),  # ***
        'HS 2012': (441700, 442110, 442110,
                    442190),  # *** - 442190 subdivided into 442191.99 in HS 2017 (442191 of bamboo)
        'HS 2017': (441700, 442110, 
                    442110,  # 442110.99 restructured for new subheading 4421.20 has been created to provideseparately for coffins in HS 2022
                    442199),
        'HS 2022': (441700, 442110, 442120, 442199)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Builder’s joinery and carpentry of wood',),
        'FAO Code': ('135',),
        'HS 1996': (441810, 441810, 441820, 441820, 441890, 441840, 441850, 441830, 441830, 441830,
                    441890, 441890, 441860, 441890, 441890),  # ***
        'HS 2002': (441810, 441810, 441820, 441820, 
                    441890,  # 441860 created based on 441890 in HS 2007
                    441840, 441850, 
                    441830, 441830, 441830,  # 441871.72.79 created based on 441830 in HS 2007
                    441890, 441890, 441860, 441890, 441890),  # ***
        'HS 2007': (441810, 441810, 441820, 441820, 441860, 441840, 441850, 441871, 441872,
                    441879, 441890, 441890, 441860, 441890, 441890),  # ***
        'HS 2012': (441810, 441810, 441820, 441820, 
                    441860,
                    441840, 441850, 
                    441871,  # new structure of 4418.7X | 441871 partially transfered to 441874 in HS 2017
                    441872,  # new structure of 4418.7X | 441872 partially transfered to 441875 in HS 2017
                    441879,  # *** - 441879 partially trasnfered into 441873 (of bamboo) in HS 2017
                    441890, 441890, 441860, 441890, 441890),  # *** - 441890 subdivided into 441891.99 in HS 2017 (441891 of bamboo)
        'HS 2017': (441810, 441810,  # 441810 subdivided into 441811.19 in HS 2022
                    441820, 441820,  # 441820 subdivided into 441821.29 in HS 2022
                    441860,  # transfer from 441860 to 441830 in HS 2022
                    441840, 441850, 441874, 441875, 441879, 
                    441899, 441899, 441860, 441899,  # transfer of products from 441860.91.99 to 4418.8X in HS 2022
                    441899),  # transfer of products from 441899 to 4418.8X and 441892
        'HS 2022': (441811, 441819, 441821, 441829, 441830, 441840, 441850, 441874, 441875, 441879, 441881, 441882, 441883, 441889, 441899)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Wooden furniture',),
        'FAO Code': ('136',),
        'HS 1996': (940130, 940140,  # ***
                    940161, 940169, 
                    940191,  # ***
                    940330, 940340, 940350, 940360, 
                    940390),  # ***
        'HS 2002': (940130, 940140,  # ***
                    940161, 940169, 
                    940190,  # ***
                    940330, 940340, 940350, 940360, 
                    940390),  # ***
        'HS 2007': (940130, 940140,  # ***
                    940161, 940169, 
                    940190,  # ***
                    940330, 940340, 940350, 940360, 
                    940390),  # ***
        'HS 2012': (940130, 940140,  # ***
                    940161, 940169, 
                    940190,  # ***
                    940330, 940340, 940350, 940360, 
                    940390),  # ***
        'HS 2017': (940130,  # *** - 940130 subdivided into 940131.39 in HS 2022
                    940140,  # *** - 940140 subdivided into 940141.49 in HS 2022
                    940161, 940169, 
                    940190,  # *** - 940190 subdivided into 940191.99 in HS 2022
                    940330, 940340, 940350, 940360, 
                    940390),  # *** - 940390 subdivided into 940391.99 in HS 2022
        'HS 2022': (940131, 940141, 940161, 940169, 940191, 940330, 940340, 940350, 940360, 940391)
    },

    {
        'FAO Product Agg': 'SECONDARY WOOD PRODUCTS',
        'FAO Code Agg': '13',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Prefabricated buildings of wood',),
        'FAO Code': ('137',),
        'HS 1996': (940600,),  # ***
        'HS 2002': (940600,),  # ***
        'HS 2007': (940600,),  # ***
        'HS 2012': (940600,),  # *** - 94.06 has been subdivided to provide separately for prefabricated buildings of wood
        'HS 2017': (940610,),
        'HS 2022': (940610,)
    },

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('',),
        'FAO Code': ('',),
        'HS 1996': (4807, 
                    481110, 481141, 481149, 481160, 481190,  # *** - 482311.19.90 partially transfered to 4811.XX in HS 2002
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482360, 482370, 482390,  # *** - complete restructuration
                    481840),  # ***
        'HS 2002': (4807, 
                    481110, 481141, 481149, 481160, 481190,  # *** - deletion of heading 4815 partially melted in 481110.90 in HS 2007 due to low volume of trade
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482360, 482370, 482390,  # *** - 482360 subdivided 482361.69 in HS 2007 (482361 of bamboo)
                    481840),  # ***
        'HS 2007': (4807, 
                    481110, 481141, 481149, 481160, 481190, 
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482369, 482370, 482390,
                    481840),  # *** - 9619 created in HS 2012 based on numerous products
        'HS 2012': (4807, 
                    481110, 481141, 481149, 481160, 481190, 
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482369, 482370, 482390,
                    9619),
        'HS 2017': (4807, 
                    481110, 481141, 481149, 481160, 481190, 
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482369, 482370, 482390,
                    9619),
        'HS 2022': (4807, 
                    481110, 481141, 481149, 481160, 481190, 
                    4814, 4816, 4817, 4818, 4819, 4820, 4821, 4822,
                    482320, 482340, 482369, 482370, 482390,
                    9619)
    },

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Composite paper and paperboard',),
        'FAO Code': ('141',),
        'HS 1996': (480710, 480790),
        'HS 2002': (480700, 480700),
        'HS 2007': (480700, 480700),
        'HS 2012': (480700, 480700),
        'HS 2017': (480700, 480700),
        'HS 2022': (480700, 480700)
    },

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Special coated paper and pulp products',),
        'FAO Code': ('142',),
        'HS 1996': (481110, 481141, 481149, 481160, 481190),  # *** - transfers of 482311.90 into 4811.XX in HS 2002 from amendments regarding the size of paper, paperboard, cellulose wadding and webs of cellulose fibres in strips or rolls
        'HS 2002': (481110, 481141, 481149, 481160, 481190),  # *** - deletion of heading 4815 because of the low volume of trade -> partially in 481110.90 in HS 2007
        'HS 2007': (481110, 481141, 481149, 481160, 481190),
        'HS 2012': (481110, 481141, 481149, 481160, 481190),
        'HS 2017': (481110, 481141, 481149, 481160, 481190),
        'HS 2022': (481110, 481141, 481149, 481160, 481190)
    },

    # N.B.: For carbon paper and copying paper ready for use, we don't take into account
    # products 481610.30 because of low volume of trade and to keep a clean classification
    # correspondence between products codes.

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Carbon paper and copying paper ready for use',),
        'FAO Code': ('143',),
        'HS 1996': (481620, 481690),  
        'HS 2002': (481620, 481690),  # 481610.30 deleted because of low volume of trade in HS 2007 (merged with 481690)
        'HS 2007': (481620, 481690),
        'HS 2012': (481620, 481690),
        'HS 2017': (481620, 481690),
        'HS 2022': (481620, 481690)
    },

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Household and sanitary paper ready for use',),
        'FAO Code': ('144',),
        'HS 1996': (481810, 481820, 481830, 481850, 481890, 
                    481840),  # ***
        'HS 2002': (481810, 481820, 481830, 481850, 481890, 
                    481840),  # ***
        'HS 2007': (481810, 481820, 481830, 481850, 481890, 
                    481840),  # *** - creation of new heading 96.19 for sanitary towels (pads) and tampons, napkins and napkin liners for babies and similar articles, of any material in HS 2012 (mix of many other codes)
        'HS 2012': (481810, 481820, 481830, 481850, 481890, 
                    961900),
        'HS 2017': (481810, 481820, 481830, 481850, 481890, 
                    961900),
        'HS 2022': (481810, 481820, 481830, 481850, 481890, 
                    961900)
    },

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Packaging cartons, boxes, etc.',),
        'FAO Code': ('145',),
        'HS 1996': (481910, 481920, 481930, 481940, 481950, 481960),
        'HS 2002': (481910, 481920, 481930, 481940, 481950, 481960),
        'HS 2007': (481910, 481920, 481930, 481940, 481950, 481960),
        'HS 2012': (481910, 481920, 481930, 481940, 481950, 481960),
        'HS 2017': (481910, 481920, 481930, 481940, 481950, 481960),
        'HS 2022': (481910, 481920, 481930, 481940, 481950, 481960)
    },

    # N.B.: For other articles of paper and paperboard ready for use, we don't take into account
    # products 481410.30 because of low volume of trade and to keep a clean classification
    # correspondence between products codes.

    {
        'FAO Product Agg': 'SECONDARY PAPER PRODUCTS',
        'FAO Code Agg': '14',
        'FAO 1982': ('n.a.',),
        'FAO Product': ('Other articles of paper and paperboard ready for use',),
        'FAO Code': ('146',),
        'HS 1996': (481420, 481490,  # ***
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482360, 482370, 482390),  # *** - lot of retsructuring on 4823.XX in HS 2002
        'HS 2002': (481420, 481490,  # *** - deletion of subheading 481430 because of low volume of trade -> merged in 481490 in HS 2007
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482360, 482370,  # *** - 482360 subdvided into 482361.69 in HS 2007 (482361 of bamboo)
                    482390),  # *** - deletion of subheading 481500 because of the low volume of trade -> merged into 482390 in HS 2007
        'HS 2007': (481420, 481490,  # *** - deletion of subheading 481410 because of low volume of trade -> merged in 481490 in HS 2012
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482369, 482370, 482390),
        'HS 2012': (481420, 481490, 
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482369, 482370, 482390),
        'HS 2017': (481420, 481490, 
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482369, 482370, 482390),
        'HS 2022': (481420, 481490, 
                    481710, 481720, 481730, 
                    482010, 482020, 482030, 482040, 482050, 482090,
                    482110, 482190, 
                    482210, 482290,
                    482320, 482340, 482369, 482370, 482390)
    }
]

correspondence_classification_pl = (
    pl.DataFrame(correspondence_classification)
        .explode('HS 2022', 'HS 2017', 'HS 2012', 'HS 2007', 'HS 2002', 'HS 1996')
        .explode('FAO Product', 'FAO Code')
        .explode('FAO 1982')
        .cast({cs.starts_with('HS'): pl.String})
)

with pl.Config(tbl_cols=-1):
    logging.info(f"\nCorrespondence table:\n{correspondence_classification_pl}")

# Saving correspondence classification
correspondence_classification_pl.write_json(
    snakemake.output[0],
)
