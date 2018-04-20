#---------------------------------------------------------------------
# Title:         IRIS - XLSX Return
# Author:        Brandon Monier
# Created:       2018-04-19 17:10:54 CDT
# Last Modified: 
#---------------------------------------------------------------------

iris_excel <- function(
    wd, geo_ser, sam_out_01, sam_out_02, sam_out_03, sam_out_04,
    sam_out_05
) {

    # Copy, timestamp, and load metadata template
    file.copy("data/seq-template.xlsx", paste0("tmp/", wd))
    setwd(paste0("tmp/", wd))

    wd2 <- substring(wd, 5)
    meta <- paste0("metadata-", wd2, ".xlsx")

    file.rename(
        from = "seq-template.xlsx", 
        to = meta
    )

    wb <- loadWorkbook(meta)


    # Metadata 01 - The first entry...
    meta_01a <- c(
        "# High-throughput sequencing metadata template (version 2.1).",
        "# All fields in this template must be completed.",
        paste(
            "# Templates containing example data are found in the METADATA",
            "EXAMPLES spreadsheet tabs at the foot of this page."
        ),
        paste(
            "# Field names (in blue on this page) should not be edited.",
            "Hover over cells containing field names to view field content",
            "guidelines."
        ),
        paste(
            "# Human data. If there are patient privacy concerns regarding",
            "making data fully public through GEO, please submit to NCBI's",
            "dbGaP (http://ww.ncbi.nlm.nih.gov/gap/) database. dbGaP has",
            "controlled access mechanisms and is an appropriate resource for",
            "hosting sensitive patient data."
        ),
        "",
        "SERIES",
        "# This section describes the overall experiment"
    )

    # Get number of contributers
    contrib <- geo_ser[grep("^geo_contrib_", names(geo_ser))]

    # Turn series into vector
    geo_ser <- unlist(geo_ser, use.names = FALSE)
    geo_ser <- c(geo_ser, rep("", 4), "title")

    # Series title info
    ser_00_title <- c(
        "title",
        "summary",
        "overall design",
        rep("contributor", length(contrib)),
        "supplementary file",
        "SRA_center_name_code"
    )

    # The second chunk of metadata info...
    meta_02a <- c(
        "",
        "SAMPLES",
        paste(
            "# This section lists and describes each of the biological", 
            "Samplesunder investgation, as well as any protocols that are",
            "specific to individual Samples."
        ),
        paste(
            "# Additional \"processed data file\" or \"raw file\" columns",
            "may be included"
        ),
        "Sample name"
    )

    # Sample `n` titles
    geo_samples <- lapply(seq_len(length(sam_out_01[[1]])), function(i) {
        paste("Sample", i)
    })
    geo_samples <- unlist(geo_samples)

    # Sample user titles
    geo_sam_title <- sam_out_01[1]
    geo_sam_title <- unlist(geo_sam_title, use.names = FALSE)
    geo_ser <- c(geo_ser, geo_sam_title)

    # Sample source name
    geo_sam_source <- sam_out_01[2]
    geo_sam_source <- unlist(geo_sam_source, use.names = FALSE)
    geo_sam_source <- c(
        rep("", length(geo_ser)),
        "source name",
        geo_sam_source
    )

    # Sample organism
    geo_sam_org <- sam_out_01[3]
    geo_sam_org <- unlist(geo_sam_org, use.names = FALSE)
    geo_sam_org <- c(
        rep("", length(geo_ser)),
        "organism",
        geo_sam_org
    )

    # Sample characteristics
    for (i in seq_len(length(sam_out_02))) {
        sam_out_02[[i]] <- c(
            rep("", length(geo_ser)),
            "characteristic",
            sam_out_02[i]
        )
    }

    # Sample molecule
    geo_sam_mol <- sam_out_03[1]
    geo_sam_mol <- unlist(geo_sam_mol, use.names = FALSE)
    geo_sam_mol <- c(
        rep("", length(geo_ser)),
        "molecule",
        geo_sam_mol
    )

    # Sample description
    geo_sam_desc <- sam_out_03[1]
    geo_sam_desc <- unlist(geo_sam_desc, use.names = FALSE)
    geo_sam_desc <- c(
        rep("", length(geo_ser)),
        "description",
        geo_sam_desc
    )

    # Sample proc. data files (main)
    geo_sam_pdf_main <- sam_out_04
    for (i in seq_len(length(geo_sam_pdf_main))) {
        geo_sam_pdf_main[i] <- paste0(geo_sam_pdf_main[i], ".txt")
    }
    geo_sam_pdf_main <- c(
        rep("", length(geo_ser)),
        "processed data file",
        geo_sam_pdf_main
    )

    # Sample proc. data files (extras)
    test_05 <- sam_out_05
    # if(all(unlist(test_05) != "")) {
        for (i in seq_len(length(sam_out_05))) {
            sam_out_05[[i]] <- c(
                rep("", length(geo_ser)),
                "processed data file",
                sam_out_05[i]
            )
        }
    # } else {
    #     NULL
    # }  




    # MASTER COLUMN 01
    master_A <- c(
        meta_01a,
        ser_00_title,
        meta_02a,
        geo_samples
    )
    # MASTER COLUMN 02
    master_B <- c(
        geo_ser
    )

    # Edit data
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = master_B,
        startRow = 9,
        startCol = 2
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_source,
        startRow = 1,
        startCol = 3
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_org,
        startRow = 1,
        startCol = 4
    )
    for (i in seq_len(length(sam_out_02))) {
        writeData(
            wb,
            sheet = "METADATA TEMPLATE ",
            x = unlist(sam_out_02[[i]], use.names = FALSE),
            startRow = 1,
            startCol = 4 + i
        )
    }
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_mol,
        startRow = 1,
        startCol = 4 + length(sam_out_02) + 1
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_desc,
        startRow = 1,
        startCol = 4 + length(sam_out_02) + 2        
    )
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = geo_sam_pdf_main,
        startRow = 1,
        startCol = 4 + length(sam_out_02) + 3      
    )
   
    if(!all(unlist(test_05, use.names = FALSE) == "")) {
        for (i in seq_len(length(sam_out_05))) {
            writeData(
                wb,
                sheet = "METADATA TEMPLATE ",
                x = unlist(sam_out_05[[i]], use.names = FALSE),
                startRow = 1,
                startCol = 4 + length(sam_out_02) + 3 + i
            )
        }
    } else {
        NULL
    }

    # Run last
    writeData(
        wb,
        sheet = "METADATA TEMPLATE ",
        x = master_A,
        startRow = 1,
        startCol = 1
    )
    # Styles...
    bodyStyle1 <- createStyle(
        fgFill = "#CCFFCC"
    )
    bodyStyle2 <- createStyle(
        wrapText = FALSE
    )
    addStyle(
        wb, 
        sheet = 1, 
        bodyStyle2, 
        rows = 1:50, 
        cols = 1:50, 
        gridExpand = TRUE
    )    
    addStyle(
        wb, 
        sheet = 1, 
        bodyStyle1, 
        rows = c(1:6, 8), 
        cols = 1:50, 
        gridExpand = TRUE
    )
    saveWorkbook(wb, meta, overwrite = TRUE)
    # return(wb)

    # setwd("../..")

    # return(paste(cat("It worked:\n\n"), test_05))
    # return(all(unlist(test_05, use.names = FALSE) != ""))
}