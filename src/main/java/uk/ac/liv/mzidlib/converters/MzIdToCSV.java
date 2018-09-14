package uk.ac.liv.mzidlib.converters;

import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzidentml.MzIdentMLParser;
import umich.ms.fileio.filetypes.mzidentml.jaxb.standard.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MzIdToCSV {

    private MzIdentMLType mzIdentMLType;
    private HashMap<String, String> fileLocationMap = new HashMap<>();
    private HashMap<String, String[]> peptideIDToSequenceMap = new HashMap<>();
    private HashMap<String, String[]> peptideRefToEvidenceMap = new HashMap<>();
    private HashMap<String, String> dbRefToName = new HashMap<>();
    private ArrayList<String> scoreName = new ArrayList<>();
    private String seq = "\t";
    private String header =  "Raw data location" + seq + "Spectrum ID"+ seq + "Spectrum Title"+ seq + "Retention Time (s)"+ seq +
            "PSM_ID"+ seq + "rank"+ seq + "Pass Threshold"+ seq + "Calc m/z"+ seq + "Exp m/z"+ seq + "Charge"+ seq + "Sequence"+ seq + "Modifications";

    public static void main(String[] args){

        try {
            new MzIdToCSV(new File(args[0]));
        } catch (FileParsingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public MzIdToCSV(File mzIdFile) throws FileParsingException, IOException {

        mzIdentMLType = MzIdentMLParser.parse(Paths.get(mzIdFile.getAbsolutePath()));

        getFileLocationMap();
        getDbRefToName();
        getPeptideIDToSequenceMap();
        getPeptideRefToEvidenceMap();
        getAllOtherPara();

        for (String name : scoreName){
            header += seq + name;
        }
        header += seq + "proteinacc_start_stop_pre_post_;"+ seq + "Is decoy";

        parseMzId(mzIdFile);

    }

    public void parseMzId(File mzIDFile) throws IOException {

        File outPutFile = new File(mzIDFile.getAbsolutePath() + ".csv");

        FileWriter fileWriter = new FileWriter(outPutFile);

        fileWriter.write(header);
        fileWriter.write("\n");

        HashMap<String, String> scoreNameToValue;

        List<SpectrumIdentificationListType> spectrumIdentificationListTypeList = mzIdentMLType.getDataCollection().getAnalysisData().getSpectrumIdentificationList();

        for (SpectrumIdentificationListType spectrumIdentificationListType : spectrumIdentificationListTypeList){
            List<SpectrumIdentificationResultType> spectrumIdentificationResultTypeList = spectrumIdentificationListType.getSpectrumIdentificationResult();

            for (SpectrumIdentificationResultType spectrumIdentificationResultType : spectrumIdentificationResultTypeList){

                String spectrumID = spectrumIdentificationResultType.getSpectrumID();

                String spectrumTitle = "";

                Double rt = -1.0;

                String spectrumFileRef = spectrumIdentificationResultType.getSpectraDataRef();

                String currentFile = fileLocationMap.get(spectrumFileRef);

                List<AbstractParamType> abstractParamTypeList = spectrumIdentificationResultType.getParamGroup();

                for (AbstractParamType abstractParamType : abstractParamTypeList){

                    if (abstractParamType instanceof CVParamType){
                        CVParamType cvParamType = (CVParamType) abstractParamType;
                        if (cvParamType.getAccession().equals("MS:1001114") || cvParamType.getAccession().equals("MS:1000016"))   {
                            rt = Double.valueOf(cvParamType.getValue());
                        }

                        if (cvParamType.getAccession().equals("MS:1000796")) {
                            spectrumTitle = cvParamType.getValue();
                        }
                    }
                }

                List<SpectrumIdentificationItemType> spectrumIdentificationItemList = spectrumIdentificationResultType.getSpectrumIdentificationItem();

                for (SpectrumIdentificationItemType spectrumIdentificationItemType : spectrumIdentificationItemList){

                    Double caMassToCharge = spectrumIdentificationItemType.getCalculatedMassToCharge();
                    Double exMassToCharge = spectrumIdentificationItemType.getExperimentalMassToCharge();
                    String psmID = spectrumIdentificationItemType.getId();
                    Integer rank = spectrumIdentificationItemType.getRank();
                    String passThresHold = String.valueOf(spectrumIdentificationItemType.isPassThreshold());
                    Integer charge = spectrumIdentificationItemType.getChargeState();
                    String peptideId = spectrumIdentificationItemType.getPeptideRef();

                    String[] peptideArray = peptideIDToSequenceMap.get(peptideId);
                    String[] evidenceArray = peptideRefToEvidenceMap.get(spectrumIdentificationItemType.getPeptideEvidenceRef().get(0).getPeptideEvidenceRef());

                    scoreNameToValue = new HashMap<>();

                    List<AbstractParamType> abstractParamTypeInList = spectrumIdentificationItemType.getParamGroup();

                    for (AbstractParamType abstractParamType : abstractParamTypeInList){
                        scoreNameToValue.put(abstractParamType.getName(), abstractParamType.getValue());
                    }

                    fileWriter.write(currentFile + "\t");
                    fileWriter.write(spectrumID + "\t");
                    fileWriter.write(spectrumTitle + "\t");
                    fileWriter.write(rt + "\t");
                    fileWriter.write(psmID + "\t");
                    fileWriter.write(rank + "\t");
                    fileWriter.write(passThresHold + "\t");
                    fileWriter.write(caMassToCharge + "\t");
                    fileWriter.write(exMassToCharge + "\t");
                    fileWriter.write(charge + "\t");
                    fileWriter.write(peptideArray[0] + "\t");
                    fileWriter.write(peptideArray[1] + "\t");
                    for (String name : scoreName){
                        if (scoreNameToValue.containsKey(name)){
                            fileWriter.write(scoreNameToValue.get(name) + "\t");
                        } else {
                            fileWriter.write("" + "\t");
                        }
                    }
                    fileWriter.write(evidenceArray[0] + "\t");
                    fileWriter.write(evidenceArray[1]);
                    fileWriter.write("\n");


                }


            }
        }
        fileWriter.close();

    }

    private void getFileLocationMap(){
        List<SpectraDataType> spectraDataTypeList = mzIdentMLType.getDataCollection().getInputs().getSpectraData();

        for (SpectraDataType spectraDataType : spectraDataTypeList){
            fileLocationMap.put(spectraDataType.getId(), spectraDataType.getLocation());
        }
    }

    private void getPeptideIDToSequenceMap(){
        List<PeptideType> peptideTypeList = mzIdentMLType.getSequenceCollection().getPeptide();

        String[] peptideArray;

        for (PeptideType peptideType : peptideTypeList){

            peptideArray = new String[2];
            peptideArray[0] = peptideType.getPeptideSequence();
            StringBuilder allMod = new StringBuilder();
            List<ModificationType> modificationTypeList = peptideType.getModification();

            for (ModificationType modificationType : modificationTypeList){
                int location = modificationType.getLocation();
                double monoMassDelta = modificationType.getMonoisotopicMassDelta();
                List<CVParamType> cvParamTypes = modificationType.getCvParam();
                String modName;

                if (cvParamTypes != null){
                    CVParamType firstType = cvParamTypes.get(0);

                    if (firstType.getName() != null){
                        modName = firstType.getName();
                    } else {
                        modName = String.valueOf(monoMassDelta);
                    }
                } else {
                    modName = String.valueOf(monoMassDelta);
                }

                allMod.append(modName).append(":").append(location);

                if (modificationTypeList.indexOf(modificationType) != modificationTypeList.size() - 1){
                    allMod.append(";");
                }
            }

            peptideArray[1] = String.valueOf(allMod);

            peptideIDToSequenceMap.put(peptideType.getId(), peptideArray);
        }

    }

    private void getDbRefToName(){

        List<DBSequenceType> dbSequenceTypeList = mzIdentMLType.getSequenceCollection().getDBSequence();

        for(DBSequenceType dbSequenceType : dbSequenceTypeList){

            dbRefToName.put(dbSequenceType.getId(), dbSequenceType.getName());

        }
    }

    private void getPeptideRefToEvidenceMap(){

        List<PeptideEvidenceType> peptideEvidenceTypeList = mzIdentMLType.getSequenceCollection().getPeptideEvidence();
        String[] peptideEvidenceArray;
        StringBuilder details;

        for (PeptideEvidenceType peptideEvidenceType : peptideEvidenceTypeList){

            peptideEvidenceArray = new String[2];
            details = new StringBuilder();

            details.append(dbRefToName.get(peptideEvidenceType.getDBSequenceRef()));
            details.append("_");
            details.append(peptideEvidenceType.getStart());
            details.append("_");
            details.append(peptideEvidenceType.getEnd());
            details.append("_");
            details.append(peptideEvidenceType.getPre());
            details.append("_");
            details.append(peptideEvidenceType.getPost());

            peptideEvidenceArray[0] = String.valueOf(details);
            peptideEvidenceArray[1] = String.valueOf(peptideEvidenceType.isIsDecoy());


            peptideRefToEvidenceMap.put(peptideEvidenceType.getId(), peptideEvidenceArray);

        }

    }

    private void getAllOtherPara(){

        for (SpectrumIdentificationListType spectrumIdentificationListType : mzIdentMLType.getDataCollection().getAnalysisData().getSpectrumIdentificationList()){
            for (SpectrumIdentificationResultType spectrumIdentificationResultType : spectrumIdentificationListType.getSpectrumIdentificationResult()){
                for (SpectrumIdentificationItemType spectrumIdentificationItemType : spectrumIdentificationResultType.getSpectrumIdentificationItem()){

                    for (AbstractParamType abstractParamType : spectrumIdentificationItemType.getParamGroup()){

                        String name = abstractParamType.getName();

                        if (!scoreName.contains(name)){
                            scoreName.add(name);
                        }

                    }
                    break;

                }
            }
        }

    }



}
