import React from 'react';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

export default function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Actions
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Show Map...</DropdownToggle>
        <DropdownMenu>
          {
            props.maps.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onMapChoice(
                  choice
                )}}
              >
                {choice.name}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
    </DropdownMenu>
  </UncontrolledDropdown>)
}