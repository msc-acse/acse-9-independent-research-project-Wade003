<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <Output_filename>
    <string_value lines="1">Test</string_value>
  </Output_filename>
  <Executable>
    <string_value lines="1">~/Packages/MultiFluids_Dev/bin/icferst</string_value>
  </Executable>
  <Input_file>
    <string_value lines="1">3d_test.mpml</string_value>
  </Input_file>
  <opal_operation>
    <importance_map/>
  </opal_operation>
  <functional>
    <end_time/>
    <location_of_interest>
      <real_value shape="3" rank="1">0.0 0.05 0.05</real_value>
    </location_of_interest>
  </functional>
  <Field_to_study name="PhaseVolumeFraction">
    <field_type>
      <scalar_field/>
    </field_type>
    <perturbations_to_do>
      <integer_value rank="0">10</integer_value>
    </perturbations_to_do>
    <Initial_condition/>
    <Gram-Schmidt/>
    <smoothing>
      <integer_value rank="0">4</integer_value>
    </smoothing>
    <Sigma>
      <real_value rank="0">0.01</real_value>
    </Sigma>
    <Phases_to_perturb>
      <integer_value shape="2" rank="1">1 2</integer_value>
    </Phases_to_perturb>
  </Field_to_study>
</OPAL>
