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
          <eigh_method/>
          <field_name name="Velocity">
            <nPOD>
              <integer_value rank="0">10</integer_value>
            </nPOD>
          </field_name>
        </svd_type>
      </compression>
      <training>
        <LSTM/>
      </training>
    </nirom>
  </opal_operation>
</OPAL>
