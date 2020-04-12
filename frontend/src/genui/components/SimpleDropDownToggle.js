import React from 'react';
import { ButtonDropdown, Col, DropdownItem, DropdownMenu, DropdownToggle, Row } from 'reactstrap';

export default function SimpleDropDownToggle(props) {
  const [dropdownOpen, setOpen] = React.useState(false);
  const [title, setTitle] = React.useState(props.title);

  return (
    <React.Fragment>
      <Row>
        <Col sm={12}>
          <props.message {...props}/>
        </Col>
      </Row>
      <Row>
        <Col sm={12}>
          <ButtonDropdown style={{width: "100%"}} isOpen={dropdownOpen} toggle={() => setOpen(!dropdownOpen)}>
            <DropdownToggle caret color="primary" block>
              {title}
            </DropdownToggle>
            <DropdownMenu style={{width: "100%"}}>
              <DropdownItem header>{props.header}</DropdownItem>
              {
                props.items.map(item => <DropdownItem onClick={() => {props.onSelect(item);setTitle(item.name)}} key={item.id}>{item.name}</DropdownItem>)
              }
            </DropdownMenu>
          </ButtonDropdown>
        </Col>
      </Row>
    </React.Fragment>
  )
}