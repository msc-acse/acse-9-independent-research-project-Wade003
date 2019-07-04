<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <Output_filename>
    <string_value lines="1">canyon</string_value>
  </Output_filename>
  <Executable>
    <string_value lines="1">~/Packages/fluidity/bin/fluidity</string_value>
  </Executable>
  <Input_file>
    <string_value lines="1">2d_canyon.flml</string_value>
  </Input_file>
  <opal_operation>
    <data_assimilation>
      <fwd_input_file>
        <string_value lines="1">2d_canyon_fwd.flml</string_value>
      </fwd_input_file>
    </data_assimilation>
  </opal_operation>
  <functional>
    <end_time/>
    <location_of_interest>
      <real_value shape="2" rank="1">1.0 0.4</real_value>
    </location_of_interest>
  </functional>
  <Time_window>
    <real_value rank="0">1</real_value>
  </Time_window>
  <random_seed/>
  <Field_to_study name="Tracer">
    <field_type>
      <scalar_field/>
    </field_type>
    <perturbations_to_do>
      <integer_value rank="0">20</integer_value>
    </perturbations_to_do>
    <Initial_condition/>
    <Gram-Schmidt/>
    <smoothing>
      <integer_value rank="0">7</integer_value>
    </smoothing>
    <use_G>
      <integer_value rank="0">1</integer_value>
    </use_G>
    <Sigma>
      <real_value rank="0">0.01</real_value>
    </Sigma>
    <Phases_to_perturb>
      <integer_value shape="1" rank="1">1</integer_value>
    </Phases_to_perturb>
  </Field_to_study>
</OPAL>
