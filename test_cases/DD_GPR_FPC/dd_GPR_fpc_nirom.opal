<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <opal_operation>
    <nirom>
      <snapshots_location>
        <path>
          <string_value lines="1">snapshots</string_value>
        </path>
        <input_file>
          <string_value lines="1">circleNRe3900_fixed.mpml</string_value>
        </input_file>
      </snapshots_location>
      <compression>
        <svd_type>
          <DD_eigh_method>
            <num_sub_base_2>
              <integer_value rank="0">2</integer_value>
            </num_sub_base_2>
          </DD_eigh_method>
          <field_name name="Velocity">
            <nPOD>
              <integer_value rank="0">10</integer_value>
            </nPOD>
          </field_name>
        </svd_type>
      </compression>
      <training>
        <DD_GPR>
          <scaling_bounds>
            <real_value shape="2" rank="1">0 10</real_value>
          </scaling_bounds>
          <constant_value>
            <real_value rank="0">1</real_value>
          </constant_value>
          <constant_bounds>
            <real_value shape="2" rank="1">1e-3 1e3</real_value>
          </constant_bounds>
          <RBF_length_scale>
            <real_value rank="0">100</real_value>
          </RBF_length_scale>
          <RBF_length_scale_bounds>
            <real_value shape="2" rank="1">1e-2 1e2</real_value>
          </RBF_length_scale_bounds>
        </DD_GPR>
      </training>
    </nirom>
  </opal_operation>
</OPAL>
