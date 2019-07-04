<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <Output_filename>
    <string_value lines="1">square</string_value>
  </Output_filename>
  <Executable>
    <string_value lines="1">/home/cheaney/Packages/MultiFluids_Dev/bin/icferst</string_value>
  </Executable>
  <Input_file>
    <string_value lines="1">het_square.mpml</string_value>
  </Input_file>
  <opal_operation>
    <importance_map/>
  </opal_operation>
  <functional>
    <end_time/>
    <location_of_interest>
      <real_value shape="2" rank="1">0.8 0.7</real_value>
    </location_of_interest>
  </functional>
  <Time_window>
    <real_value rank="0">0.1</real_value>
  </Time_window>
  <Field_to_study name="PhaseVolumeFraction">
    <field_type>
      <scalar_field/>
    </field_type>
    <perturbations_to_do>
      <integer_value rank="0">50</integer_value>
    </perturbations_to_do>
    <Initial_condition/>
    <Gram-Schmidt/>
    <smoothing>
      <integer_value rank="0">5</integer_value>
    </smoothing>
    <use_G>
      <integer_value rank="0">4</integer_value>
    </use_G>
    <Sigma>
      <real_value rank="0">0.1</real_value>
    </Sigma>
    <Phases_to_perturb>
      <integer_value shape="1" rank="1">1</integer_value>
    </Phases_to_perturb>
  </Field_to_study>
</OPAL>
